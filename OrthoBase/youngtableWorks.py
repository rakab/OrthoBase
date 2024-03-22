import z3

from .utils import *
database = {}

class YoungTables(object):
    def __init__(self, obj):
        self.tables = obj

    def __mul__(self,other):
        print("Mult ",self,other)
        results = []
        for tab in self.tables:
            results.append(tab*other)
        lst = []
        for res in results:
            for tab in res.tables:
                lst.append(tab)
        return(self.__class__(lst))

class YoungTable(object):
    Nc = 4
    def __init__(self, obj):
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
            self.__init__(part_rows)
            #self.dims = (nrows,ncols)
            #self.table = obj
            #lst=[obj[i] for i in obj if obj[i]!=0]
            #if lst:
                #self.ncells = len(lst)
            #else:
                #self.ncells = 0
            ##for col in ncols:
                ##if row in
        self.dim = self.calc_dim()

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
        for row in range(self.dims[0]):
            for col in range(self.dims[1]):
                #if self.table[col,row] != 0:
                print(self.table[row,col], end=' ')
            print()
        print("Mult:", self.dim)

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


    def __mul__(self,other):
        if self.dim < other.dim:
            return(other.__mul__(self))
        sol = database.get((self.dim,other.dim))
        if sol:
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
        print(rows,cols)
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

            solutions.append(self.__class__(table))
            getDifferentSolutionMatrix(sol,mod,y,rows,cols)
            #break
        print("# Solutions: ", num_solutions)
        database[(self.dim,other.dim)] = solutions
        return(YoungTables(solutions))
