#! /usr/bin/python
"""
__author__: Iva Pritisanac

MAGMA engine - executes methods defined in classes of magma_class.py, and graph_heuristics.py;
"""
import sys,os
from magma_class import Magma
from graph_heuristics import Graph

# import the user input file with parameters for the calculation
"""
try:
    filename = sys.argv[1]
except:
    print "ERROR: expecting input file name as the first command line argument"
    sys.exit(2)

if len(sys.argv)==3:
    magma_version = sys.argv[2]
else:
    magma_version = "c"    # default will be c version

if magma_version=="c" or magma_version=="py":
    print "Running MAGMA version %s"%magma_version
else:
    print "ERROR in the input MAGMA version"
    print "Set to the default: c version"
    magma_version="c"
"""
   
magma_version="py"
filename = "input_files/hsp90/input_hsp90.txt"
    
M = Magma(filename) # make an instance of Magma class
run_mode,minsize,mcesmode,stripmode,niter,nitermces,runtime,check_isomorphism = M.parse_variables(filename) # run parser to get variables
labels,data_vertices,data_adj,struct_vertices,struct_adj,ligand=M.parse_data(filename) # run parser to get data structures

G = Graph(labels,labels,'type',ligand)  # on the basis of the input parameters instantiate Graph for subsequent methods
data_vertices,data_adj,struct_vertices,struct_adj=M.check_graphs(data_vertices,data_adj,struct_vertices,struct_adj,G,minsize)

if run_mode=='all_subgraphs':
    print "Running MAGMA including all data subgraphs"
    if check_isomorphism and M.subgraph_isomorphism(data_vertices,data_adj,struct_vertices,struct_adj,G,"_all"):
        print "Calculation finished. Clean program exit ..."
        sys.exit(0)
    else:       
        print "Subgraph isomorphism not found \nRunning MCES algorithm..."
        runtag = "full"
        assign_options=M.set_assignment_options(data_vertices,data_adj,struct_vertices,struct_adj,G,stripmode,runtag,filter_dict={})
        # update initial data structures according to the results of the optimization runs
        assign_options,data_vertices,data_adj,runtime= M.optimise_run_order(G,assign_options,data_vertices,data_adj,struct_vertices,struct_adj,niter,nitermces,runtime,runtag,version=magma_version,optimise_mode="y")
        if magma_version=="c":  # if c version of mcgregor algorithm is employed
            M.run_complete_mcgregor_c(assign_options,data_vertices,data_adj,struct_vertices,struct_adj,runtime,runtag,mcesmode)
        if magma_version=="py": # if python version of mcgregor algorithm is employed
            M.run_complete_mcgregor_py(assign_options,data_vertices,data_adj,struct_vertices,struct_adj,runtime,runtag,mcesmode)
    
    print "Calculation finished. Clean program exit ..."
    M.analyse_magma_results()
    sys.exit(0)
        
elif run_mode=='connected_subgraphs':
    
    print "Running MAGMA in the split subgraphs mode"
    non_isomorphic=[]
    non_isomorphic_id=[]

    subgraph_tag=0
    connected_subgraphs = M.split_data_graph(data_vertices,data_adj,G)    # split NMR data graph to its connected subgraphs
    for subgraph in connected_subgraphs:
        subgraph_tag+=1
        subgraph_vert,subgraph_adj = M.get_subgraph_data(subgraph,G)
        if check_isomorphism and M.subgraph_isomorphism(subgraph_vert,subgraph_adj,struct_vertices,struct_adj,G,str(subgraph_tag)):
            print "Assigned ", str(subgraph_tag), " with VF2 algorithm"
        else:
            print "At short distance threshold, ", str(subgraph_tag)," data subgraph is not isomorphic to any structure subgraph"
            non_isomorphic.append(subgraph) # collect subgraphs that could not be evaulated with VF2 for mces evaluation
            non_isomorphic_id.append(subgraph_tag)  # collect ids of subgraphs that could not be evaluated
    
    if len(non_isomorphic)>0:   # if there are any subgraphs left that need MCES evaluation
        print "Non subgraph isomorphic = ",len(non_isomorphic)
        for stag,left_subgraph in enumerate(non_isomorphic):
            runtag = str(stag)
            sub_vert,sub_adj = M.get_subgraph_data(left_subgraph,G)
            # //TO DO// support filtering here
            #filter_assign_options = M.get_filter(subgraph_vertices,subgraph_adj,struct_vertices,struct_adj,stripmode)
            #set_assignment_options(...filter=filter_assign_options)
            sub_assign_options=M.set_assignment_options(sub_vert,sub_adj,struct_vertices,struct_adj,G,stripmode,runtag,filter_dict={})
            # update initial data structures according to the results of the optimization runs
            sub_assign_options,sub_vert,sub_adj,runtime = M.optimise_run_order(G,sub_assign_options,sub_vert,sub_adj,struct_vertices,struct_adj,niter,nitermces,runtime,runtag,version=magma_version,optimise_mode="y")
            if magma_version=="c":
                M.run_complete_mcgregor_c(sub_assign_options,sub_vert,sub_adj,struct_vertices,struct_adj,runtime,runtag,mcesmode)
            if magma_version=="py":
                M.run_complete_mcgregor_py(sub_assign_options,sub_vert,sub_adj,struct_vertices,struct_adj,runtime,runtag,mcesmode)
    
    print "Exiting magma split!"
    print "Calculation finished. Clean program exit ..."
    M.analyse_magma_results()
    sys.exit(0)
else:
    print "ERROR: unknown run mode found in the input file %s\t%s"%(run_mode,filename)
    print "Exiting program ..."
    sys.exit(1)
