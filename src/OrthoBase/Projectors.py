import os
import sys
import shutil
from importlib.resources import files as imp_files

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

    def run(self):
        f_new = open(os.path.join(self.path,"new_projs.frm"), "w")
        f_new.write(imp_files('OrthoBase.FORMData').joinpath('header.frm').read_text())
        f_new.write('\n')
        for m in self.mults:
            #print("Constructing projector for ", m.dim)
            m.decompose()
            #print("First occurence",m.first_occ)
            #n=m
            #gen=1
            #while n.parent1 is not None:
            #    gen+=1
            #    n=n.parent1
            #print("Generation",gen)
            if self.mults.gen==m.first_occ:
                f_new.write(self.proj_new(m))
        f_new.close()

    def yng_proj(self, y, i):
        """
        Construct a non-hermitian Young operator
        """
        print("Projecting ",y.dim,y.part_cols,y.part_rows,y.table)
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
        y1 = m.decomposition[0]
        y2 = m.decomposition[1]
        print("New proj {0}={1}x{2}".format(m.dim,m.parent1.dim,8),y1.dim,y2.dim)
        m.print()
        proj_name = "P{0}(".format(m.dim)
        oldproj1 = "P{0}(".format(m.parent1.dim)
        colfac1 = ""
        for g in range(1,self.mults.gen):
            proj_name += "mu{0},".format(g)
            proj_name += "nu{0},".format(g)
            oldproj1 += "mu{0},".format(g)
            oldproj1 += "mmu{0},".format(g)
            colfac1 += "SUNT(mmu{0},i{0},j{0})*".format(g)
        proj_name += "mu{0},".format(self.mults.gen)
        proj_name += "nu{0},".format(self.mults.gen)
        proj_name = proj_name[:-1]
        proj_name += ")"
        oldproj1 = oldproj1[:-1]
        oldproj1 += ")"
        oldproj2 = oldproj1.replace("mu","nu")
        colfac1 += "SUNT(mu{0},i{0},j{0})".format(self.mults.gen)
        colfac2 = colfac1.replace("mu","nu")
        colfac2 = colfac2.replace("i","ri")
        colfac2 = colfac2.replace("j","rj")
        colfac2 = colfac2.replace("ri","jj")
        colfac2 = colfac2.replace("rj","ii")
        y1_code=self.yng_proj(y2, "i")
        y2_code=self.yng_proj(y2, "j")
        #Subtract the rest
        subtr = "("
        for pM in m.decomposition[2]:
            subtr += "P{0}+".format(pM.dim)
        subtr = subtr[:-1]
        subtr += ")"
        proj = "L {0}={1}*{2}*{3}*{4}*{5}*{6}-{7};\n".format(proj_name,oldproj1,colfac1,oldproj2,colfac2,y1_code,y2_code,subtr)
        print(proj)
        return(proj)
        #m.decomposition[1].print()
        #print(m.first_occ)
