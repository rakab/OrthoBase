import os
import sys
import shutil
import subprocess
from itertools import product
from importlib.resources import files as imp_files
from mpi4py import MPI
from mpi4py.futures import MPIPoolExecutor
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from .utils import *

class Projectors(object):
    def __init__(self, multiplets, path):
        self.mults = multiplets
        if os.path.exists(path):
            print("{0} folder exists...".format(path))
            if os.access(path, os.W_OK | os.X_OK):
                if any(os.scandir(path)):
                    user_ans=input("Directory is not empty. The content will be overwritten. Proceed (y/[n])? ")
                    if user_ans.lower() == "y":
                        for root, dirs, files in os.walk(path):
                            for f in files:
                                os.unlink(os.path.join(root, f))
                            for d in dirs:
                                shutil.rmtree(os.path.join(root, d))
                    else:
                        print("Exiting...")
                        sys.exit(0)
                self.path = path
            else:
                raise PermissionError("No write access for {0} directory!".format(path))
        else:
            raise FileNotFoundError("{0} directory does not exist!".format(path))
        #Copy SUn.prc to the output dir
        shutil.copy(imp_files('OrthoBase.FORMData').joinpath('SUn.prc'),self.path)
        print("########")

    @property
    def nodes(self):
        return(self._nodes)

    @nodes.setter
    def nodes(self, value):
        if not isinstance(value,list):
            raise TypeError("nodes must be a list of strings!")
        self._nodes = value

    @property
    def parallel_evaluation(self):
        return(self._parallel_evaluation)

    @parallel_evaluation.setter
    def parallel_evaluation(self,value):
        if not isinstance(value, bool):
            raise TypeError("parallel_evaluation must hold a boolean value!")
        self._parallel_evaluation = value

    @property
    def FORM_np(self):
        return(self._FORM_np)

    @FORM_np.setter
    def FORM_np(self,value):
        if not isinstance(value, int):
            raise TypeError("FORM_np must be an integer number of processes!")
        self._FORM_np = value

    def run(self):
        logger.info("Building expressions of projectors")
        f_old = open(os.path.join(self.path,"old_projs.frm"), "w")
        f_old_list = open(os.path.join(self.path,"old_projs_list.frm"), "w")
        f_new = open(os.path.join(self.path,"new_projs.frm"), "w")
        f_new_list = open(os.path.join(self.path,"new_projs_list.frm"), "w")
        f_new.write(imp_files('OrthoBase.FORMData').joinpath('header.frm').read_text())
        f_new.write('\n')
        f_old.write(imp_files('OrthoBase.FORMData').joinpath('header.frm').read_text())
        f_old.write('\n')
        oldproj_names = "CF "
        clebsch_names = "CF "
        clebsches = set()
        for m in self.mults:
            m.decompose()
            if self.mults.gen==m.first_occ:
                proj, oldproj = self.proj_new(m)
                f_new_list.write(proj)
                oldproj_names += oldproj + ","
            else:
                clebsch,rhs_proj_name,projs = self.proj_old(m)
                f_old_list.write("\n".join(projs))
                f_old_list.write("****************\n")
                clebsches.add(clebsch)
                clebsches.add(rhs_proj_name)
        clebsch_names += ",".join(clebsches)
        f_old.write(clebsch_names+";\n\n")
        f_old.write("#include old_projs_list.frm\n\n")
        f_new.write(oldproj_names[:-1]+";\n\n")
        f_new.write("#include new_projs_list.frm\n\n")
        #f_new.write("L test = SUNT;")
        f_new.write(imp_files('OrthoBase.FORMData').joinpath('footer_new.frm').read_text())
        f_old.write(imp_files('OrthoBase.FORMData').joinpath('footer_old.frm').read_text())
        f_new_list.close()
        f_new.close()
        f_old.close()
        f_old_list.close()

        logger.info("Starting to run FORM")
        os.chdir(self.path)
        process=subprocess.Popen(["tform",f"-w{self.FORM_np}","new_projs.frm"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
        print("Program output:")
        stdout_lines = []
        while True:
            line = process.stdout.readline()
            if not line:
                break
            stdout_lines.append(line)
        process.wait()

        stdout = "".join(stdout_lines)

        print(stdout)


    def yng_operator(self, y, i):
        """
        Construct a non-hermitian Young operator
        """
        #print("Projecting ",y.dim,y.part_cols,y.part_rows,y.table)
        logger.debug(f"Projecting {y.dim},{y.part_cols},{y.part_rows},{y.table}")
        code = ""
        #symmetric combinations
        #lhs contains the indices going inside the diagram
        #rhs contains the indices coming out of the diagram
        for row,cols in enumerate(y.part_cols):
            lhs = "perm_(f1,"
            rhs = "f2("
            idx = sum(y.part_cols[:row])
            for col in range(cols):
                idx += 1
                #symm += "d(j{0},mj{0}),".format(idx)
                lhs += "{0}{1},".format(i,idx)
                rhs += "m{0}{1},".format(i,idx)
            lhs = lhs[:-1]
            rhs = rhs[:-1]
            lhs += ")"
            rhs += ")"
            code += "{0}*{1}*".format(lhs,rhs)
        #asymmetric combinations
        for col,rows in enumerate(y.part_rows):
            lhs = "perm_(1,f1,"
            rhs = "f2("
            for row in range(rows):
                idx = sum(y.part_cols[:row])+col+1
                #asymm += "d(mj{0},jj{0}),".format(idx)
                lhs += "m{0}{1},".format(i,idx)
                rhs += "{0}{0}{1},".format(i,idx)
            lhs = lhs[:-1]
            rhs = rhs[:-1]
            lhs += ")"
            rhs += ")"
            code += "{0}*{1}*".format(lhs,rhs)
        code = code[:-1]
        return code


    def proj_new(self, m):
        """
        Construction of the projectors which appear for the first time in the
        creation history

        Notation:
        i indices - quarks
        j indices - antiquarks
                             Y
                           /  \
                          /i1  \ii1
        mu1 -- P -- mmu1--       ------ mnu1 -- P -- nu1
         .                \j1   /jj1                  .
         .                 \   /                      .
         .                   Yb                       .
         .                                            .
         .               /       \                    .
        mu.... ---------   ......  ----------------- nu...
                         \       /
        """
        y_qb = m.decomposition[0]
        y_q = m.decomposition[1]
        if y_qb.ncells != m.first_occ:
            y_qb = y_qb.conjugate()
        elif y_q.ncells != m.first_occ:
            y_q = y_q.conjugate()
        print("New proj {0}={1}x{2} decomposed as {3}x{4}".format(m.dim_txt,m.parent1.dim_txt,8,y_qb.dim,y_q.dim))
        m.print()
        #parents=",".join(obj.dim_txt for obj in reversed((lambda x: [x] + (lambda x: [x] + [x.parent1])(x.parent1) if x.parent1 else [x])(m)) if obj)
        parents=",".join(obj.dim_txt for obj in reversed(m.parent_list()))+","+m.dim_txt
        proj_name = "[P({0})]".format(parents)
        all_indices = "("
        #oldproj1 = "P[{0}](".format(m.parent1.dim)
        oldproj_names = "CF "
        oldproj = "[P({0})]".format(parents[:parents.rindex(',')])
        oldproj1 = oldproj + "("
        colfac1 = ""
        for g in range(1,self.mults.gen):
            all_indices += "mu{0},".format(g)
            all_indices += "nu{0},".format(g)
            oldproj1 += "mu{0},".format(g)
            oldproj1 += "mmu{0},".format(g)
            colfac1 += "SUNT(mmu{0},i{0},j{0})*".format(g)
        all_indices += "mu{0},".format(self.mults.gen)
        all_indices += "nu{0},".format(self.mults.gen)
        all_indices = all_indices[:-1]
        all_indices += ")"
        oldproj1 = oldproj1[:-1]
        oldproj1 += ")"
        oldproj2 = oldproj1.replace("mu","nu")
        colfac1 += "SUNT(mu{0},i{0},j{0})".format(self.mults.gen)
        colfac2 = colfac1.replace("mu","nu")
        colfac2 = colfac2.replace("i","ri")
        colfac2 = colfac2.replace("j","rj")
        colfac2 = colfac2.replace("ri","jj")
        colfac2 = colfac2.replace("rj","ii")
        y_q_code=self.yng_operator(y_q, "i")
        y_qb_code=self.yng_operator(y_qb, "j")
        #Subtract the rest
        subtr = "("
        for pM in m.decomposition[2]:
            subtr += "[P({0})]".format(parents[:parents.rindex(',')]+","+str(pM.dim))
            oldproj += ",[P({0})]".format(parents[:parents.rindex(',')]+","+str(pM.dim))
            subtr += f"{all_indices}+"
        subtr = subtr[:-1]
        subtr += ")"
        proj_name += all_indices
        proj = "L {0}={1}*{2}*{3}*{4}*{5}*{6}-{7};\n".format(proj_name,oldproj1,colfac1,oldproj2,colfac2,y_q_code,y_qb_code,subtr)
        print(proj)
        return(proj,oldproj)
        #m.decomposition[1].print()
        #print(m.first_occ)

    def proj_old(self, m):
        """
        Construction of the projectors which previously already appeared in the
        creation history
        """

        print("Old proj {0}={1}x{2}".format(m.dim_txt,m.parent1.dim_txt,8))
        parents=",".join(obj.dim_txt for obj in reversed(m.parent_list()))+","+m.dim_txt
        proj_name = "[C({0})]".format(parents)
        rhs_proj_name = "[P({0})]".format(m.dim)
        clebsch = "[C({0})]".format(parents[:parents.rindex(',')])
        clebsch_idx = "("

        for g in range(1,self.mults.gen):
            clebsch_idx += "mu{0},".format(g)

        colfac1 = ""
        m.parent1.decompose()
        i_list = list()
        j_list = list()
        for g in range(1,m.parent1.first_occ+1):
            clebsch_idx += "mmu{0},".format(g)
            colfac1 += "SUNT(mmu{0},i{0},j{0})*".format(g)
            i_list.append(f"i{g}")
            j_list.append(f"j{g}")
        colfac1 += "SUNT(mu{0},i{0},j{0})".format(self.mults.gen)
        i_list.append(f"i{self.mults.gen}")
        j_list.append(f"j{self.mults.gen}")
        clebsch_idx = clebsch_idx[:-1]
        clebsch_idx += ")"

        ii_list = list()
        jj_list = list()
        colfac2 = ""
        rhs_idx = "("
        if m.first_occ==0:
            #colfac2 = "1"
            colfac2 = "d_(nu1,nu2)"
            rhs_proj_name = "1"
            rhs_idx = ""
        else:
            for g in range(1,m.first_occ+1):
                rhs_idx += "mnu{0},".format(g)
                rhs_idx += "nu{0},".format(g)
                colfac2 += "SUNT(mnu{0},jj{0},ii{0})*".format(g)
                ii_list.append(f"ii{g}")
                jj_list.append(f"jj{g}")
            colfac2 = colfac2[:-1]
            rhs_idx = rhs_idx[:-1]+")"

        perms_ij  = list(product(i_list,j_list))
        perms_iii = list(product(i_list,ii_list))
        perms_iijj= list(product(ii_list,jj_list))
        perms_jjj = list(product(j_list,jj_list))
        all_perms = perms_ij+perms_iii+perms_iijj+perms_jjj
        logger.debug(f"Indices to pair: i_list={i_list}, j_list={j_list}, ii_list={ii_list}, jj_list={jj_list}")
        pairs=generate_pairings(all_perms,len(i_list)+len(ii_list))
        results = list()
        for n,pairing in enumerate(pairs):
            proj_name_n = proj_name.replace("]", f"{n}]")
            connector='*'.join(f'd_({i},{j})' for i, j in pairing)
            proj = f"L {proj_name_n}={clebsch}{clebsch_idx}*{colfac1}*{colfac2}*{connector}*{rhs_proj_name}{rhs_idx};\n"
            results.append(proj)
            print(proj)
        if rhs_proj_name == "1":
            rhs_proj_name = ""
        return(clebsch,rhs_proj_name,results)
