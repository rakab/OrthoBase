import z3
#import pickle
import copy
import re
from collections import Counter
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from .utils import *
database = {}

class YoungTableaux(object):
    def __init__(self, obj, g=1):
        self.tables = obj
        dims=set()
        for sol in self.tables:
            self.gen = g+1
            dims.add(sol.dim_txt)
        self.dims = dims

    def __getitem__(self, i):
        if isinstance(i,int):
            return self.tables[i]
        elif isinstance(i,str):
            return [tab for tab in self.tables if tab.dim_txt == i]

    def __len__(self):
        return len(self.tables)

    def __mul__(self,other):
        print("Mult ",self,other)
        results = []
        for tab in self.tables:
            results.append(tab*other)
        lst = []
        for res in results:
            for tab in res.tables:
                lst.append(tab)
        #YoungTableaux.gene += 1
        return(self.__class__(lst,self.gen))

    def print(self,latex=False):
        print(self.dims)
        for d in self.dims:
            mults={tab for tab in self.tables if tab.dim_txt == d}
            unique_shapes = len(set(tuple(m.part_rows) for m in mults))
            if unique_shapes > 1:
                logger.info(f"Renaming ambiguous multiplets..")
                for m in mults:
                    #print(m.part_cols)
                    m.dim_txt += "_"+f"({','.join(str(x) for x in m.part_cols)})"
        counts = Counter([sol.dim_txt for sol in self.tables])
        pattern = r'(\d+)(?:b)?(?:_\(\d+(?:,\d+)*\))?'
        #sorted_tables = sorted(counts.items(), key=lambda item: ((lambda s: int(s[:-1]) if s.endswith('b') else int(s))(item[0])))
        logger.info(f"Sorting multiplets for printing...")
        sorted_tables = sorted(counts.items(), key=lambda item: ((lambda s: int(re.search(pattern, s).group(1)))(item[0])))
        print("Total number of tables:", len(self.tables))
        print("Decomposition:")
        if latex:
            latex_code=""
            normal=""
            for mult, count in sorted_tables:
                normal += "{0}*{1}+".format(count,mult)

                mult_name = re.search(r'(\d+)', mult).group(1)
                has_b = 'b' in mult
                subscript = re.search(r'_\((.*?)\)', mult)

                # Build the LaTeX string
                if has_b:
                    latex_code += f'{count}\\cdot\\overline{{{mult_name}}}'
                else:
                    latex_code += f'{count}\\cdot {mult_name}'

                if subscript:
                    latex_code += f'_{{({subscript.group(1)})}}'
                latex_code += r" \oplus "
            latex_code = (lambda s: s.rsplit(' ', 1)[0])(latex_code[:-1])
            print(normal[:-1])
            print("Latex source:")
            print(latex_code)
        else:
            normal=""
            for mult, count in sorted_tables:
                normal += "{0}*{1}+".format(count,mult)
            print(normal[:-1])

    def __iter__(self):
        return iter(self.tables)

class YoungTableau(object):
    @staticmethod
    def init_from_list(obj,Nc):
        while obj[-1] == 0:
            obj.remove(0)
        while obj[0] == Nc and len(obj)>1:
            obj.remove(Nc)
        part_rows = obj
        dims = (max(part_rows),len(part_rows))
        table = {}
        ncells = 0
        for row in range(dims[0]):
            for col in range(dims[1]):
                if part_rows[col] > row:
                    table[(row,col)] = 1
                    ncells += 1
                else:
                    table[(row,col)] = 0
        part_cols = list()
        for i in range(part_rows[0]):
            part_cols.append(sum(1 for j in part_rows if j>i))
        return (part_rows,part_cols,dims,ncells,table)

    def __init__(self, obj, Nc, parent1=None, parent2=None):
        self.Nc = Nc
        self.init_parents(parent1,parent2)
        #part_cols - (partition by colons) how many colons are in each row
        #part_rows - (partition by rows) how many rows are in each column
        if isinstance(obj,list):
            self.part_rows,self.part_cols,self.dims,self.ncells,self.table = self.init_from_list(obj,self.Nc)
        elif isinstance(obj,dict):
            nrows = list(obj)[-1][0]+1
            ncols = list(obj)[-1][1]+1
            part_rows = list()
            for col in range(ncols):
                n = 0
                for row in range(nrows):
                    if obj[(row,col)] != 0:
                        n = n+1
                part_rows.append(n)
            self.__init__(part_rows,self.Nc,parent1,parent2)
        self.dim = (self.calc_dim())
        conj_rows = list()
        for r in self.part_rows:
            conj_rows.append(self.Nc-r)
        conj_rows.append(self.Nc)
        conj_rows.reverse()
        conj_part_cols = self.init_from_list(conj_rows,self.Nc)[1]
        if len(self.part_cols)>len(conj_part_cols):
            self.dim_txt = str(self.dim)+"b"
        elif len(self.part_cols)==len(conj_part_cols) and self.part_cols[-1] > conj_part_cols[-1]:
            self.dim_txt = str(self.dim)+"b"
        else:
            self.dim_txt = str(self.dim)

    def init_parents(self,parent1,parent2):
        if parent1 is not None and parent1.parent1 is not None:
            self.parent1 = parent1
            self.parent2 = parent2
        else:
            self.parent1 = parent2
            self.parent2 = parent1

    def calc_dim(self):
        """
        Calculate the dimension of the multiplet
        """
        num = 1
        den = 1
        rows = self.dims[0]
        cols = self.dims[1]
        for i in range(rows):
          for j in range(cols):
            val = self.table[(i,j)]
            hooklen = 0
            for jj in range(j+1,cols):
              if self.table[(i,jj)] != 0: hooklen+=1
            for ii in range(i+1,rows):
              if self.table[(ii,j)] != 0: hooklen+=1
            if val  != 0:
              #print(val, end=' ')
              num *= self.Nc+j-i
              den *= hooklen+1
        return(int(num/den))

    def print(self):
        #print(self.part_cols)
        print(f"Y={self.part_rows} Dimension: {self.dim_txt}")
        for row in range(self.dims[0]):
            for col in range(self.dims[1]):
                #if self.table[col,row] != 0:
                print(self.table[row,col], end=' ')
            print()
        for col in range(self.dims[1]):
            print(u'\u2500', end=' ')
        print()

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index >= len(self.table):
            raise StopIteration
        else:
            item = list(self.table)[self._index]
            position = self._index+1
            self._index += 1
            if self.table[item] == 0:
                raise StopIteration
            return item, position

    def relabel(self):
        table = list()
        for row in range(self.dims[0]):
            for col in range(self.dims[1]):
                if self.table[row,col] != 0:
                    table.append(row*self.dims[1]+col+1)
        table.append(0)
        return(table)

    def id2labels(self,x):
        return((x-1)/self.dims[1]+1)

    def conjugate(self):
        conj_rows = list()
        for r in self.part_rows:
            conj_rows.append(self.Nc-r)
        conj_rows.append(self.Nc)
        conj_rows.reverse()
        return(self.__class__(conj_rows,self.Nc))

    def parent_list(self):
        parents = list()
        p = self.parent1
        while p != None:
            parents.append(p)
            p = p.parent1
        return(parents)

    def decompose(self):
        """
        Decomposition of a table in terms of quark and anti-quark projection
        operators as described in Appendix B of arXiv:1207.0609 + some
        optimization
        """
        logger.debug(f"Decomposing of {self.dim}-{self.part_rows} into q qbar tables")
        conj_rows = list()
        for r in self.part_rows:
            conj_rows.append(self.Nc-r)
        l = int(len(self.part_rows)/2)
        direction = 0 #Helps to avoid looping back and forth
        m=0
        n=0
        while l>0 and l<len(self.part_rows):
            m = sum(conj_rows[:l])
            n = sum(self.part_rows[l:])
            if m < n:
                new_direction = 1
                l += 1
            elif m > n:
                new_direction = -1
                l -= 1
            else:
                break
            if direction != 0 and direction != new_direction:
                break
            direction = new_direction

        if n==0:
            a = None
            b = None
            c = None
            logger.debug(f"Decomposition ended, firts occ: {n}")
        elif m!=n:
            a = None
            b = None
            c = None
            n = None
            print("This multiplet never appears in the product of adjoint representations.")
        else:
            l = len(self.part_rows)
            a = self.__class__(self.part_rows[:l-n],self.Nc)
            b = self.__class__(self.part_rows[l-n:],self.Nc)
            c=[m for m in a*b if m.dim != self.dim]
            logger.debug(f"Decomposition ended, {self.dim}={a.dim}x{b.dim} firts occ: {n}")
        self.first_occ = n
        #a and b are the antiquark and quark tables
        #c contains all the projectors which have to be subtracted
        self.decomposition = (a,b,c)

    def __copy__(self):
        copy_ = type(self).__new__(self.__class__)
        copy_.__dict__.update(self.__dict__)
        copy_.parent1 = None
        copy_.parent2 = None
        return copy_

    def __mul__(self,other):
        if self.dim < other.dim:
            return(other.__mul__(self))
        sol = database.get((tuple(self.part_cols),tuple(other.part_cols)))
        if sol:
            #sol = pickle.loads(pickle.dumps(sol,-1))
            #sol = copy.deepcopy(sol)
            for s in sol:
                s = copy.copy(s)
                s.init_parents(self,other)
            return(YoungTableaux(sol))
        block_num=99
        x = {}
        y = {}
        sol = z3.Solver()
        other_idx = other.relabel()
        rows = self.dims[0]+other.ncells
        cols = self.dims[1]+other.part_cols[0]
        if rows > self.Nc:
            rows = self.Nc
        #print("multing",self.dim,other.dim)
        #print(rows,cols)
        for row in range(rows):
            for col in range(cols):
                x[(row, col)] = z3.Int("x(%i,%i)" % (row, col))
                y[(row, col)] = z3.Int("y(%i,%i)" % (row, col))
                if row in range(self.dims[0]) and col in range(self.dims[1]):
                    if self.table[(row,col)] != 0:
                        sol.add(x[(row,col)] == block_num)
                    else:
                        #sol.add(z3.Or(x[(col, row)] == 1, x[(col, row)] == 0 ))
                        member_of(sol,x[(row,col)],other_idx)
                else:
                        #sol.add(z3.Or(x[(col, row)] == 1, x[(col, row)] == 0 ))
                        member_of(sol,x[(row,col)],other_idx)
                if col>0:
                    sol.add(z3.If(z3.And(x[(row, col-1)]==0, x[(row,col)]!=0),False,True))
                    sol.add(z3.If(z3.And(x[row,col-1]!=block_num,x[row,col-1]!=0,x[row,col]!=block_num,x[row,col]!=0),z3.If(x[row,col]>x[row,col-1],True,False),True))
                if row>0:
                    sol.add(z3.If(z3.And(x[(row-1,col)]==0, x[(row,col)]!=0),False,True))

        for row in range(rows):
            for col in range(cols):
                for row_i in range(rows):
                    for col_i in range(cols):
                        if row_i != row or col_i != col:
                            sol.add(z3.If(z3.And(x[row,col]!=block_num,x[row,col]!=0),x[(row, col)] != x[(row_i, col_i)],True))
                #Rewrite id of each cell using a,b notation (a=1,b=2,c=3,...)
                sol.add(z3.If(z3.And(x[row,col]!=block_num,x[row,col]!=0),y[(row, col)] == other.id2labels(x[(row,col)]),y[(row, col)] == x[(row, col)]))

        y_flat = [y[(row,col)] for row in range(rows) for col in reversed(range(cols))]
        repetitions = [z3.Int("p%i" % i) for i in range(other.part_rows[0])]
        #Count a,b label
        for iy,yy in enumerate(y_flat):
            for lab in range(other.part_rows[0]):
                repetitions[lab]=z3.Sum([z3.If(y_flat[i] == lab+1, 1,0) for i in range(iy+1)])
            #Reject If more b than a
            for lab1 in range(other.part_rows[0]):
                for lab2 in range(lab1,other.part_rows[0]):
                    if lab1 == lab2: continue
                    sol.add(repetitions[lab1]>=repetitions[lab2])
        #No repeating labels (a,b,c,...) in the column
        for lab in range(other.part_rows[0]):
            for col in range(cols):
                sol.add(z3.Sum([z3.If(y[(row,col)] == lab+1, 1,0) for row in range(rows)])<=1)

        for col in range(cols):
            for row in range(1,rows):
                sol.add(z3.If(z3.And(y[(row-1,col)]!=block_num,x[row,col]!=0,y[(row,col)]<y[(row-1,col)]),False,True))

        #Total num of cells = ncells in left table + ncells in right table
        sol.add(z3.Sum([z3.If(x[row,col] != 0, 1,0) for row in range(rows) for col in range(cols)]) == self.ncells + other.ncells )

        #Remove cols with zeroes
        #for row in range(rows):
            #z3.If(z3.Sum([x[(col,row)] for col in range (cols)])==0,x[(2,3)]==1,x[(2,3)]==1)

        #print(sol)
        z3.set_param("parallel.enable", True)
        #z3.set_param('parallel.threads.max', 4)
        num_solutions = 0
        solutions = list()
        while sol.check() == z3.sat:
            num_solutions += 1
            mod = sol.model()
            table = {}
            for row in range(rows):
                for col in range(cols):
                    table[(row, col)] = mod.eval(y[(row,col)]).as_long()
                    if table[(row,col)] == block_num:
                        table[(row,col)] = 9
            #remove_cols = 0
            #for row in range(rows):
                #lst=[table[(col, row)] for col in range(cols) if table[(col,row)]==0]
                #if not lst:
                    #remove_cols = row
#
            #ttable = {}
            #for col in range(cols):
                #for row in range(remove_cols+1,rows):
                    #ttable[(col,row-remove_cols-1)] = table[(col,row)]

            solutions.append(self.__class__(table,self.Nc,self,other))
            getDifferentSolutionMatrix(sol,mod,y,rows,cols)
            #break
        logger.debug(f"# Solutions: {num_solutions}")
        database[(tuple(self.part_cols),tuple(other.part_cols))] = solutions
        return(YoungTableaux(solutions))
