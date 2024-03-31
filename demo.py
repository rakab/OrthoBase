#!/usr/bin/env python3
import pickle
from OrthoBase import YoungTools as YT
from OrthoBase import Projectors as P
import logging
logging.basicConfig(level=logging.WARNING)
logging.getLogger('OrthoBase').setLevel(logging.INFO)

#Define number of colors
Nc = 3
#Define SU(3) octet (gluon)
g = YT.YoungTableau([Nc-1,1],Nc)
a = YT.YoungTableau([2,2,2,2,2,1,1,1,1,1,1,1,1],Nc)
b = YT.YoungTableau([2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1],Nc)
aa = YT.YoungTableau([2,2,2,2,2,2,2,2,1,1,1,1,1],Nc)
bb = YT.YoungTableau([2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1],Nc)
#a.print()
#b.print()

multiplets = YT.YoungTableaux([a,b,a,a,b,aa,bb])

#Perform the decomposition of gggggggggg state
#multiplets = g*g*g*g*g*g*g*g*g*g*g
#multiplets = g*g*g*g*g*g*g*g*g*g*g
#multiplets = g*g*g
#print(multiplets[1].parent1.parent1, multiplets[1].parent2.dim)
#Print information about the resulting multiplets
print("#####################")
#for m in multiplets:
    #m.print()
multiplets.print(1)
#print(multiplets["27"])

#Construct projectors
#projectors = P(multiplets, '/tmp/test/')
#projectors.parallel_evaluation = True
#projectors.nodes = [
    #"MySuperComputer1.edu",
    #"MySuperComputer2.edu",
    #"MySuperComputer3.edu",
    #]
#Number of parallel processes for MPI
#projectors.mpi_np = 900 #For simple processes it is pointless to have more processes than number of multiplets
#Number of threads FORM is allowed to use
#projectors.FORM_np = 12
#projectors.run()
