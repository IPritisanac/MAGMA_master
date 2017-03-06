#! /usr/bin/python

"""
MAGMA engine - executes methods defined in classes of modules: data, graph_heuristics, molecule, etc. (c) Iva Pritisanac | iva.pritisanac@gmail.com
"""
import sys,os
from magma_class import Magma
from mces import McGregor
from graph_heuristics import Graph
# import from user input file with parameters set for the calculation

"""
try:
    filename = sys.argv[1]
except:
    print "ERROR: expecting input file name as the first command line argument"
    sys.exit(2)

if len(sys.argv)==2:
    magma_version = sys.argv[2]
else:
    magma_version = "c"    # default will be c version
"""

filename = "input_files/input_hsp90.txt"
magma_version = "c"
   
if magma_version=="c" or magma_version=="py":
    print "Running MAGMA version %s"%magma_version
else:
    print "ERROR in the input MAGMA version"
    print "Set to the default: c version"
    magma_version="c"
    
M = Magma(filename) # make an instance of Magma class
run_mode,minsize,mcesmode,stripmode,niter,nitermces,runtime = M.parse_variables(filename) # run parser to get variables
labels,data_vertices,data_adj,struct_vertices,struct_adj,ligand=M.parse_data(filename) # run parser to get data structures

G = Graph(labels,labels,'type',ligand)  # on the basis of the input parameters instantiate Graph for subsequent methods
data_vertices,data_adj,struct_vertices,struct_adj=M.check_graphs(data_vertices,data_adj,struct_vertices,struct_adj,G,minsize)

if run_mode=='all_subgraphs':
    print "Running MAGMA including all data subgraphs"
    if M.subgraph_isomorphism(data_vertices,data_adj,struct_vertices,struct_adj,"_all",G):
        print "Calculation finished. Clean program exit ..."
        sys.exit(0)
    else:       
        print "Subgraph isomorphism not found \nRunning MCES algorithm..."
        runtag = "full"
        if magma_version=="c":  # if c version of mcgregor algorithm is employed
            assign_options=M.set_assignment_options(data_vertices,data_adj,struct_vertices,struct_adj,G,stripmode,filter_dict={})
            # update initial data structures according to the results of the optimization runs
            assign_options,data_vertices,data_adj,runtime= M.optimise_run_order(G,assign_options,data_vertices,data_adj,struct_vertices,struct_adj,niter,nitermces,runtime,runtag,version=magma_version,optimise_mode="y")
            M.run_complete_mcgregor(assign_options,data_vertices,data_adj,struct_vertices,struct_adj,runtime,runtag,mcesmode)
        if magma_version=="py": # if pythonv version of mcgregor algorithm is employed
            pass
        
    print "Calculation finished. Clean program exit ..."
    sys.exit(0)
    
elif run_mode=='connected_subgraphs':
    print "Running MAGMA in the split subgraphs mode"
    #M.magma_split()
    pass # run split magma (method of class magma)
    print "Calculation finished. Clean program exit ..."
    sys.exit(0)
else:
    print "ERROR: unknown run mode found in the input file %s\t%s"%(run_mode,filename)
    print "Exiting program ..."
    sys.exit(1)