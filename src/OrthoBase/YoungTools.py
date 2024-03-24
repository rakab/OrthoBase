import z3
#import pickle
import copy
import numpy as np
#from collections import Counter
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from .utils import *
database = {}

class YoungTables(object):
    def __init__(self, obj, g=1):
        self.tables = obj
        dims=list()
        for sol in self.tables:
            dims.append(sol.dim)
            self.gen = g+1
        #Recalculate dimensions if n and n-bar exist
        #for i,d in enumerate(dims):
            #if i == 0:
                #continue
            #for ii,dd in enumerate(dims[:i]):
                #if d != dd:
                    #continue
                ##n and n-bar are the same, e.g. octet
                #if self.tables[i].part_cols == self.tables[ii].part_cols:
                    #continue
                ##n and n-bar are different e.g. 10 and 10-bar
                #sol_tomod = None
                #if len(self.tables[i].part_cols)>len(self.tables[ii].part_cols):
                    #sol_tomod = self.tables[i]
                #elif len(self.tables[i].part_cols)<len(self.tables[ii].part_cols):
                    #sol_tomod = self.tables[ii]
                #else:
                    #if self.tables[i].part_cols[-1] > self.tables[ii].part_cols[-1]:
                        #sol_tomod = self.tables[i]
                    #else:
                        #sol_tomod = self.tables[ii]
                #if sol_tomod is not None and sol_tomod.dim[-1]!='b':
                    #sol_tomod.dim = sol_tomod.dim+'b'
        ############################################

    def __getitem__(self, i):
        return self.tables[i]
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
        #YoungTables.gene += 1
        return(self.__class__(lst,self.gen))

    def print(self):
        dims=list()
        for sol in self.tables:
            dims.append(sol.dim)
        mults, counts = np.unique(dims, return_counts=True)
        latex=""
        normal=""
        for mult, count in np.c_[mults,counts]:
            normal += "{1}*{0}+".format(mult,count)
            latex += "{1}\\cdot {0}\\oplus".format(mult,count)
        print("Total number of tables:", len(self.tables))
        print("Decomposition:")
        print(normal[:-1])

    def __iter__(self):
        return iter(self.tables)

class YoungTable(object):
    def __init__(self, obj, Nc, parent1=None, parent2=None):
        self.Nc = Nc
        self.init_parents(parent1,parent2)
        #part_cols - (partition by colons) how many colons are in each row
        #part_rows - (partition by rows) how many rows are in each column
        if isinstance(obj,list):
            while obj[-1] == 0:
                obj.remove(0)
            while obj[0] == self.Nc and len(obj)>1:
                obj.remove(self.Nc)
            self.part_rows = obj
            self.dims = (max(self.part_rows),len(self.part_rows))
            table = {}
            self.ncells = 0
            for row in range(self.dims[0]):
                for col in range(self.dims[1]):
                    if self.part_rows[col] > row:
                        table[(row,col)] = 1
                        self.ncells += 1
                    else:
                        table[(row,col)] = 0
            self.table = table
            self.part_cols = list()
            for i in range(self.part_rows[0]):
                self.part_cols.append(sum(1 for j in self.part_rows if j>i))
        elif isinstance(obj,dict):
            #for col,row in obj:
                #print(col,row)
                #if prev_col != col:
                    #print(col,prev_row)
                #prev_col = col
                #prev_row = row
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
            #self.dims = (nrows,ncols)
            #self.table = obj
            #lst=[obj[i] for i in obj if obj[i]!=0]
            #if lst:
                #self.ncells = len(lst)
            #else:
                #self.ncells = 0
            ##for col in ncols:
                ##if row in
        self.dim = (self.calc_dim())

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
        print("Dimension:", self.dim)
        for row in range(self.dims[0]):
            for col in range(self.dims[1]):
                #if self.table[col,row] != 0:
                print(self.table[row,col], end=' ')
            print()
        for col in range(self.dims[1]):
            print(u'\u2500', end=' ')
        print()

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
        l = len(self.part_rows)
        while True:
            m = sum(conj_rows[:int(l/2)])
            n = sum(self.part_rows[int(l/2):])
            if m < n:
                l = l+1
            elif m > n:
                l = l-1
            else:
                break
        if n==0:
            a = None
            b = None
            c = None
            logger.debug(f"Decomposition ended, firts occ: {n}")
        else:
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
        for attr in self.__dict__:
            copy_.__dict__[attr] = self.__dict__[attr]
            #if attr in ["parent1","parent2"]:
                # Make an actual copy of the attribute
                #pass
                #copy_.__dict__[attr] = copy.copy(self.__dict__[attr])
            #else:
                # Copy reference
                #copy_.__dict__[attr] = self.__dict__[attr]
        return copy_

    def __deepcopy__(self,memo):
        #return YoungTable(self.table,self.Nc,self.parent1,self.parent2)
        copy_ = type(self).__new__(self.__class__)
        for attr in self.__dict__:
            if attr in ["parent1","parent2"]:
                # Make an actual copy of the attribute
                copy_.__dict__[attr] = copy.copy(self.__dict__[attr])
            else:
                # Copy reference
                copy_.__dict__[attr] = self.__dict__[attr]
        return copy_

    def __mul__(self,other):
        if self.dim < other.dim:
            return(other.__mul__(self))
        sol = database.get((tuple(self.part_cols),tuple(other.part_cols)))
        if sol:
            #sol = pickle.loads(pickle.dumps(sol,-1))
            sol = copy.deepcopy(sol)
            for s in sol:
                s.init_parents(self,other)
            return(YoungTables(sol))
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
        return(YoungTables(solutions))
