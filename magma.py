#! /usr/bin/python
# MAGMA engine (c) author:Iva Pritisanac | iva.pritisanac@gmail.com
# executes methods defined in mces.py

import sys,os,copy
from mces import McGregor	# import all methods of the classes in mces.py
from data import ParseInput,PDBData,NMRData
from graph_heuristics import Graph,Heuristic
from subgraph_isomorphism import IgraphSubIso

import numpy as np
import igraph as ig
import networkx as nx
import time
import random

#filename = sys.argv[1]	# input text file
filename = "input_hsp90.txt"
# create instances of all classes

P=ParseInput()
try:
    P.parse(filename)
except Exception, e:
    print "Could not parse input file:\n%s"%e
    sys.exit(0)
try:
    P.check_variables()
except Exception, e:
    print "ERROR in input variables!\n"
    print e
    sys.exit(0)
#    Read in PDB file and generate appropriate graph (as specified by the input parameters)
try:
    PDB = PDBData(P.pdb_file,P.chains,P.residues,P.distance_type,P.short_distance_threshold,P.long_distance_threshold,P.include_ligand,P.ligand_chain,P.ligand_name,P.ligand_atoms)
except Exception, e:
    print "Error when creating instance of a class PDBData"
    print e
    sys.exit(0)
#    Read in NMR files and generate appropriate graph (as specified by the input parameters)
try:
    NMR = NMRData(P.HMQC_list,P.format_peak_list,P.labelled_residues,P.read_NOEs_from_file,P.NOE_list,P.format_NOE_list,P.dimensions_NOEs,P.read_extracted_NOEs,P.extracted_NOE_list,P.protein_ligand_noe_file,P.check_NOE_reciprocity,P.check_NOE_intensity,P.merge_NOEs,P.merge_LV_label,P.include_ligand_noes)

except Exception, e:
    print "Error when creating instance of a class NMRData"
    print e
    sys.exit(0)
###########################################################################################################################
positions,extended_info,info = PDB.ExtractPositionsIndices(P.pdb_flag)
dist_table_short,info_short,dist_matrix_short,mask_dist_short = PDB.InteractionNetwork(positions,info,P.short_distance_threshold,P.dist_matrix,P.indx_file)
dist_table_long,info_long,dist_matrix_long,mask_dist_long = PDB.InteractionNetwork(positions,info,P.long_distance_threshold,P.dist_matrix,P.indx_file)
if P.include_ligand:
    print "Adding ligand connections to the protein structure graph!"
    lig_dist_table,lig_info,lig_dist_matrix,lig_dist_mask=PDB.LigandInteractionNetwork(positions,extended_info,P.ligand_distance_threshold)
pdb_dict_short = PDB.GetConnectionDict(dist_table_short,info_short,False)
pdb_dict_long = PDB.GetConnectionDict(dist_table_long,info_long,False)
#PDB.WriteDistanceFile(pdb_dict_long_distances,"dist.txt")
#PDB.PrettyPlotContacts(dist_table_short, positions, dist_matrix_short, mask_dist_short, P.short_distance_threshold)
#PDB.PrettyPlotContacts(dist_table_long, positions, dist_matrix_long, mask_dist_long, P.long_distance_threshold)
###########################################################################################################################
## if we do not discriminate Leucines from Valines and vice versa
if P.merge_LV_label:
    print "Leucine and Valine labels are not known\nMerging..."
    pdb_dict_short = PDB.TurnValToLeu(pdb_dict_short)
    pdb_dict_long = PDB.TurnValToLeu(pdb_dict_long)
else:
    print "Leucines and Valines labels are known!"
if P.include_ligand:
    pdb_dict_lig = PDB.GetConnectionDict(lig_dist_table,lig_info,True)
    pdb_dict_short= PDB.JoinLigandConnections(pdb_dict_short,pdb_dict_lig)
    pdb_dict_long= PDB.JoinLigandConnections(pdb_dict_long,pdb_dict_lig)
    short_node_list,short_adjacency = PDB.GetNodesAdjacency(pdb_dict_short,True)
    long_node_list,long_adjacency = PDB.GetNodesAdjacency(pdb_dict_long,True)
else:
    short_node_list,short_adjacency = PDB.GetNodesAdjacency(pdb_dict_short)
    long_node_list,long_adjacency = PDB.GetNodesAdjacency(pdb_dict_long)
############################################################################################################################
NMR.CheckInputFormats()
NMR.InputDimensions()
NMR.SetNOEDict()
if P.include_ligand_noes:
    noe_dict = NMR.JoinLigandConnections(NMR.noes,NMR.lig_noes)
    noe_node_list,noe_adjacency = NMR.GetNodesAdjacency(noe_dict,True)
else:
    noe_node_list,noe_adjacency = NMR.GetNodesAdjacency(NMR.noes)    
##############################################################################################################################
if P.include_ligand_noes:
    G = Graph(P.labelled_residues,P.labelled_residues,"type",True)
else:
    G = Graph(P.labelled_residues,P.labelled_residues,"type")
##    Remove disconnected nodes from the G1 - smaller graph
noe_node_list,noe_adjacency = G.RemoveDisconnectedNodes(noe_node_list,noe_adjacency)
noe_node_list,noe_adjacency = G.RemoveSmallSubgraphs(noe_node_list,noe_adjacency,P.min_size)
############################################################################################################################################
# Basic graphs checks
#
small_graph = G.NetworkxGraph(noe_node_list,noe_adjacency)
big_graph = G.NetworkxGraph(short_node_list,short_adjacency)
#############################################################################################################################################
#G.PlotGraphLabels(small_graph,"eiso1.pdf")
G.CheckGraphsSize(small_graph,big_graph)
G.CheckGraphsLabels(small_graph,big_graph)
print "Basic NOE graph analysis:\n"
print "N (nodes, NOE graph) = ",len(small_graph.nodes())
print "N (edges, NOE graph) = ",len(small_graph.edges())
print "N (nodes, PDB graph) = ", len(big_graph.nodes())
print "N (edge, PDB graph) = ", len(big_graph.edges())
print "NOE sparsity >> ", len(small_graph.edges())/float(len(big_graph.edges()))
print "average degree >> ",np.mean(nx.degree(small_graph).values())
print "median degree >>",np.median(nx.degree(small_graph).values())
#############################################################################################################################################
# Splits NOE connections to separate files for separate connected components
start_time = time.time()
run_mces_joined = False    # flag that indicates if MCES algorithm will be run on all subgraphs together
n_conn_subgraphs = 0    # counter will be identifier of connected components of the NOE graph
subgraph_isomorphic = 0 # counter of subgraph isomorphic connected components
t_subgraph_isomorphic = time.time()
subgraph_non_isomorphic = 0 # counter of non subgraph isomorphic components of NOE graph
other_subgraphs = []
###############################################################################################################################################
if P.mode == "all_subgraphs":	#len(noe_node_list) < 50: # if the size of the NOE graph is sufficiently small -- check subgraph iso for all together
    graph_noe,graph_indexing = G.IgraphGraph(noe_node_list,noe_adjacency)
    graph_pdb_short,graph_pdb_short_indexing = G.IgraphGraph(short_node_list,short_adjacency)
    EP = IgraphSubIso() #    instance of a class for subgraph isomorphism check and extraction
    # test subgraph isomorphism to PDB graph at a distance threshold (short threshold)
    start_all_vf2 = time.time()
    if EP.IgraphSubIsomorphism(graph_noe,graph_pdb_short):
        n_solutions = EP.IgraphCountSubIsomorphism(graph_noe,graph_pdb_short)
        all_results = EP.IgraphListSubIsomorphism(graph_noe,graph_pdb_short,graph_indexing,graph_pdb_short_indexing,[],[],noe_node_list,[],short_node_list,[],[],False)
        labels_result_dict=EP.IgraphSubIsomorphismWriteResult(P.outfile_vf2+"complete",str(n_conn_subgraphs),all_results,graph_indexing,graph_pdb_short_indexing)
        print "SI condition met for the entire NOE graph!\nN (explained NOEs) = ",len(small_graph.edges())
        print "Took >> ", time.time()-start_all_vf2, "seconds"
        print "Finished assigning, exiting the program ... "	
        sys.exit(0)

    else:
        run_mces_joined = True # set the complete MCES run flag on!         
        candidates_dict = G.Compare(noe_node_list,noe_adjacency,short_node_list,short_adjacency,"Jaccard",False)
        H = Heuristic()
        # execute graph matching order heuristics based on the Jaccard coefficients for direct neighbourhoods and the Hungarian algorithm 
        munkres_candidates = H.HungarianFirstOrder(noe_node_list,noe_adjacency,short_node_list,short_adjacency,candidates_dict,1)
        candidates_dict = G.CombineMunkresCandidates(candidates_dict,munkres_candidates,False) 
        node_scores = {}
        node_lists = {}
        node_adjacencies = {}
        tot_n_iter = 0
        time_start = time.time()
        ##
        # to find the best starting point in NOE graph for the graph traversal
        # start at every node and build a traversal based on BFS
        # run a single iteration of MCES algorithm
        # collect a score given this choice    
        print "Starting a routine for re-ordering the vertex matching priorities"   
        while tot_n_iter<P.prior_iter:
            tot_n_iter+=1
            for n in range(len(noe_node_list)):
                start_time = time.time()
                shuffled_noe_node_list = G.OptimalGraphOrderFromNode(noe_node_list[n],noe_node_list,noe_adjacency)
                shuffled_noe_adjacency = G.ReOrderAdjecancy(shuffled_noe_node_list,noe_node_list,noe_adjacency)
                node_lists.setdefault(copy.deepcopy(noe_node_list[n]),copy.deepcopy(shuffled_noe_node_list))
                node_adjacencies.setdefault(copy.deepcopy(noe_node_list[n]),copy.deepcopy(shuffled_noe_adjacency))
                MC=McGregor(shuffled_noe_node_list)
                MC.PrepareMcGregor(shuffled_noe_adjacency,short_node_list,short_adjacency,long_node_list,long_adjacency,candidates_dict)
                # runall=True/False,n_mces=None/integer,time_check=False/True,maximum_time=None/integer
                iter_score,results = MC.McGregorLoop(P.outfile_McGregor+"_complete_"+str(n)+".txt",False,1,False,None,False)
                node_scores.setdefault(iter_score,[]).append(copy.deepcopy(noe_node_list[n]))                
            best_start_node = random.choice(node_scores[max(node_scores.keys())])	# choose randomly from a pool of the best starting vertices
            print "Started the mapping search from >> ",best_start_node
            MC=McGregor(node_lists[best_start_node])
            MC.PrepareMcGregor(node_adjacencies[best_start_node],short_node_list,short_adjacency,long_node_list,long_adjacency,candidates_dict)
            iter_score,best_results = MC.McGregorLoop(P.outfile_McGregor+"_"+"_best"+".txt",True,P.n_prior_mces,True,P.max_run_time,False)
            print "In iteration ",tot_n_iter,"scored >> ", iter_score
            candidates_dict = MC.RePrioritize(best_results,candidates_dict)
        MC=McGregor(node_lists[best_start_node]) #  an instance of McGregor class
        MC.PrepareMcGregor(node_adjacencies[best_start_node],short_node_list,short_adjacency,long_node_list,long_adjacency,candidates_dict)
        
        if P.mces_mode == "one":
            print "\nstarting exact search for a MCES...\n"            
            iter_score,best_results = MC.McGregorLoop(P.outfile_McGregor+"_"+"_final_best_one"+".txt",False,None,False,None,False)
            time_took = time.time() - time_start
            print 'Took >> ', time_took, " to find 1 MCES"
            
        elif P.mces_mode == "all":
            print "\nstarting exact search for all MCES...\n"
            iter_score,best_results = MC.McGregorLoop(P.outfile_McGregor+"_final_best_all"+".txt",True,None,False,None,False)
            time_took = time.time() - time_start
            print 'Took >> ', time_took, " to find all MCES"
   
    if run_mces_joined:
        sys.exit(0)
################################################################################################################################################
# separate over connected subgraphs
connected_subgraphs = G.SplitConnectionsOverSubgraphs(noe_node_list,noe_adjacency)
print "N connected NOE graph components >> ", len(connected_subgraphs)
# separate graph to subgraph components
for subgraph in connected_subgraphs:
    sub_conn_dict = G.GetConnectionDict(subgraph)
    if P.include_ligand and P.include_ligand_noes:
        sub_noe_node_list,sub_noe_adjacency= G.GetNodesAdjacency(sub_conn_dict,True)   
    else:
        sub_noe_node_list,sub_noe_adjacency= G.GetNodesAdjacency(sub_conn_dict)
    t_subgraph_isomorphic = time.time() # puts time limit to isomorphism testing
    ################################################################################################################################
    """
    Call implementation of Vf2 algorithm (c++) implemented in python package igraph
    """
    if P.check_subgraph_iso:    
        n_conn_subgraphs+=1
        subgraph_noe,graph_noe_indexing = G.IgraphGraph(sub_noe_node_list,sub_noe_adjacency)
        graph_pdb_short,graph_pdb_short_indexing = G.IgraphGraph(short_node_list,short_adjacency)
        start_one_vf2 = time.time()
        EP = IgraphSubIso() #    instance of a class for subgraph isomorphism check and extraction
        # for each of connected data subgraphs (of size larger than given in the input file (> ~specified N nodes))
        # test subgraph isomorphism to PDB graph at a distance threshold (short threshold)
            
        if EP.IgraphSubIsomorphism(subgraph_noe,graph_pdb_short):
            end_one_vf2 = time.time()-start_one_vf2
            print "Took >> ",end_one_vf2, " to evaluate subgraph isomorphism!"
            print "\nSubgraph isomorphism detected at short distance threshold\n"
            subgraph_isomorphism = True
            subgraph_isomorphic+=1
            n_solutions = EP.IgraphCountSubIsomorphism(subgraph_noe,graph_pdb_short)
            start_all_vf2=time.time()
            all_results = EP.IgraphListSubIsomorphism(subgraph_noe,graph_pdb_short,graph_noe_indexing,graph_pdb_short_indexing,[],[],sub_noe_node_list,[],short_node_list,[],[],False)
            end_all_vf2 = time.time()-start_all_vf2
            print "Took >>",end_all_vf2," to evaluate all SI!"
            labels_result_dict=EP.IgraphSubIsomorphismWriteResult(P.outfile_vf2,str(n_conn_subgraphs),all_results,graph_noe_indexing,graph_pdb_short_indexing)
            print "N (subgraph isomorphisms) = ", n_solutions
            print "Result written to a file >> ",P.outfile_vf2+"_"+str(n_conn_subgraphs)+".txt"
            print "Assigned",n_conn_subgraphs,"NOE subgraph!\n" 
            #other_subgraphs.append(subgraph)    # if testing for comparison mcgregor/vf2
        else:
            print "\n\nAt ",str(P.short_distance_threshold),"A >> ",n_conn_subgraphs,"NOE subgraph is not isomorphic to any pdb subgraph\n\n"
            subgraph_isomorphism = False
            subgraph_non_isomorphic+=1                
            other_subgraphs.append(subgraph)
    else:
        print "\n\nNOTE: Check for subgraph isomorphism is SWITCHED OFF; Proceeding with MCES algorithm ...\n\n"
        subgraph_isomorphism = False
        subgraph_non_isomorphic+=1                
        other_subgraphs.append(subgraph)
##################################################################################################################################################
print "Out of", len(connected_subgraphs),"connected components: "
print "N(subgraph isomorphic) = ",subgraph_isomorphic
print "N(non subgraph isomorphic) = ",subgraph_non_isomorphic
time_start = time.time()
mcsubgraph = 0
if len(other_subgraphs) > 0:
    print "Maximal common subgraph algorithm applies to ", len(other_subgraphs)," subgraph(s) ... \n"
    # for every non subgraph isomorphic subgraph
    for subgraph in other_subgraphs:
        mcsubgraph+=1
        """
        Apply heuristic before MCES search
        """
        sub_conn_dict = G.GetConnectionDict(subgraph)
        if P.include_ligand_noes:
            noe_node_list,noe_adjacency= G.GetNodesAdjacency(sub_conn_dict,True)
        else:
            noe_node_list,noe_adjacency= G.GetNodesAdjacency(sub_conn_dict)
        candidates_dict = G.Compare(noe_node_list,noe_adjacency,short_node_list,short_adjacency,"Jaccard",False)
        H = Heuristic()
        ##
        #    heuristics based on Jaccard coefficients for direct neighbourhoods and Hungarian algorithm
        munkres_candidates = H.HungarianFirstOrder(noe_node_list,noe_adjacency,short_node_list,short_adjacency,candidates_dict,1)        
        candidates_dict = G.CombineMunkresCandidates(candidates_dict,munkres_candidates,False) 	
        node_scores = {}
        tot_n_iter = 0
        ##
        # to find the best starting point in NOE graph for the graph traversal
        # start at every node and build a traversal based on BFS
        # run a single iteration of MCES algorithm
        # collect a score given this choice               
        while tot_n_iter<P.prior_iter:
            tot_n_iter+=1
            for n in range(len(noe_node_list)):
                start_time = time.time()
                shuffled_noe_node_list = G.OptimalGraphOrderFromNode(noe_node_list[n],noe_node_list,noe_adjacency)
                shuffled_noe_adjacency = G.ReOrderAdjecancy(shuffled_noe_node_list,noe_node_list,noe_adjacency)
                MC=McGregor(shuffled_noe_node_list)
                MC.PrepareMcGregor(shuffled_noe_adjacency,short_node_list,short_adjacency,long_node_list,long_adjacency,candidates_dict)
                # runall=True/False,n_mces=None/integer,time_check=False/True,maximum_time=None/integer
                iter_score,results = MC.McGregorLoop(P.outfile_McGregor+"_"+str(mcsubgraph)+"_"+str(n)+".txt",False,1,False,None,False)
                node_scores.setdefault(iter_score,[]).append(noe_node_list[n])
            
            # based on the score for all node start options, choose the best one
            best_start_node = random.choice(node_scores[max(node_scores.keys())])   # choose randomly from a pool of the best ones
            print "Started the mapping search from >> ",best_start_node
            shuffled_noe_node_list = G.OptimalGraphOrderFromNode(best_start_node,noe_node_list,noe_adjacency)   # reorder node list according to the BFS with root being the best starting node
            shuffled_noe_adjacency = G.ReOrderAdjecancy(shuffled_noe_node_list,noe_node_list,noe_adjacency) # reorder NOE adjacency according to the NOE list as given above
            # write the results of the initial search to the output file
            MC=McGregor(shuffled_noe_node_list)
            MC.PrepareMcGregor(shuffled_noe_adjacency,short_node_list,short_adjacency,long_node_list,long_adjacency,candidates_dict)
            iter_score,best_results = MC.McGregorLoop(P.outfile_McGregor+"_"+str(mcsubgraph)+"_best"+".txt",True,P.n_prior_mces,True,P.max_run_time,False)
            print "In iteration ",tot_n_iter,"scored >> ", iter_score
            candidates_dict = MC.RePrioritize(best_results,candidates_dict)
            
        shuffled_noe_node_list = G.OptimalGraphOrderFromNode(best_start_node,noe_node_list,noe_adjacency)   # final NOE nodes list
        shuffled_noe_adjacency = G.ReOrderAdjecancy(shuffled_noe_node_list,noe_node_list,noe_adjacency) #    final adjacency list
        MC=McGregor(shuffled_noe_node_list) #    instance of McGregor class
        MC.PrepareMcGregor(shuffled_noe_adjacency,short_node_list,short_adjacency,long_node_list,long_adjacency,candidates_dict)
        if P.mces_mode == "one":
            print "\nstarting exact search for a MCES...\n"            
            iter_score,best_results = MC.McGregorLoop(P.outfile_McGregor+"_"+str(mcsubgraph)+"_final_best_one"+".txt",False,None,False,None,True)
            time_took = time.time() - time_start
            print 'Took >> ', time_took, "seconds to find 1 MCES"
            
        elif P.mces_mode == "all":
            print "\nstarting exact search for all MCES...\n"
            iter_score,best_results = MC.McGregorLoop(P.outfile_McGregor+"_"+str(mcsubgraph)+"_final_best_all"+".txt",True,None,False,None,False)
            time_took = time.time() - time_start
            print 'Took >> ', time_took, "seconds to find all MCES"
        else:
            print "No mces_mode defined\n Exiting..."
            sys.exit(100)
            
"""
//To do//
Analyse solutions
Plot solutions
"""