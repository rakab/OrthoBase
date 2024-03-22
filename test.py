from YoungTools import YoungTable
import numpy as np

#a = YoungTable([2,2,1,1])
#a = YoungTable([2,2,2,2,2,1,1,1,1,1,])
#a = YoungTable([2,2,2,2,1])
#a = YoungTable([2,1,1,1,1])
#a = YoungTable([1,1,1,1,1,1])
#a = YoungTable([1,1,1,1])
#a = YoungTable([2,2,2,2,1])
#a.print()
#aa = YoungTable([2])

Nc=3
b = YoungTable([Nc-1,1],Nc)
#b = YoungTable([2])
b.print()

#bb = YoungTable([2,2,2,2,1,1])
#cc = YoungTable([2,2,2,1,1])

a = YoungTable([1],Nc)
z = YoungTable([2],Nc)
#c = b*a
#c = b*b
c = b*b
#c = b*b*b*b*b*b
#for sol in c.tables:
    #sol.print()
dims=list()
for sol in c.tables:
   #print(sol.dim,end=' ')
   #sol.print()
   dims.append(sol.dim)
print()
print(sum(dims))
#print(dims.count(1)+dims.count(8)+dims.count(10)+dims.count(27))
mults, counts = np.unique(dims, return_counts=True)
print(mults)
print("Sols: ",len(c.tables))

latex=""
for mult, count in np.c_[mults,counts]:
    latex += "{1}\\cdot {0}\\oplus".format(mult,count)

#print(latex)

print("#############################")
