#!/usr/bin/env python3
from OrthoBase import YoungTools as YT
from OrthoBase import Projectors as P
import logging
logging.basicConfig(level=logging.WARNING)
logging.getLogger('OrthoBase').setLevel(logging.INFO)

#Define the number of colors
Nc = 3
#Define the SU(3) octet (gluon)
g = YT.YoungTableau([Nc-1,1],Nc)
multiplets = g*g

#Print information about the resulting multiplets
print("#####################")
multiplets.print()

#Construct projectors
projectors = P(multiplets, '/tmp/test/')
#projectors.parallel_evaluation = True
#projectors.nodes = [
    #"MySuperComputer1.edu",
    #"MySuperComputer2.edu",
    #"MySuperComputer3.edu",
    #]
#Number of parallel processes for MPI
#projectors.mpi_np = 900 #For simple processes it is pointless to have more processes than number of multiplets
#Number of threads FORM is allowed to use
projectors.FORM_np = 12
projectors.run()
