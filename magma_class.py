#! /usr/bin/python

"""
__author__: Iva Pritisanac

executes methods defined in mces.py
"""

import sys, os, copy
from mces import McGregor	# import all methods of the class
from mces import MCES_PY
from data import ParseInput,PDBData,NMRData
from graph_heuristics import Heuristic
from subgraph_isomorphism import IgraphSubIso

import analysis

import numpy as np
import networkx as nx
import time
import random

class Magma():
    
    def __init__(self,inputf):
        
        self.dirpath=os.path.dirname(os.path.abspath(inputf)) # extract directory name out of the input filename
        if not os.path.exists(self.dirpath+os.sep+"results"): # if the output results directory does not exist on this path
            os.makedirs(self.dirpath+os.sep+"results")   # create the directory
        self.outdir=self.dirpath+os.sep+"results"   # set the output directory
        
        #sys.stdout = open(self.outdir+os.sep+'log','a')   # set screen output to a file
        
    ## run parser and return only variables -- different modes, integer parameter values etc.    
    def parse_variables(self,input_file):
        
        P=ParseInput()  # create instance of the parser class to parse the input file
        try:
            P.parse(input_file)
        except Exception, e:
            print "ERROR: Could not parse input file:\n%s"%e
            sys.exit(1)
        try:
            P.check_variables()
        except Exception, e:
            print "ERROR in input variables!\n"
            print e
            sys.exit(1)
            
        return P.mode,P.min_size,P.mces_mode,P.strip_mode,P.prior_iter,P.n_prior_mces,P.max_run_time,P.check_subgraph_iso
    
    ## run parser to fill out and return data structures
    def parse_data(self,input_file):

        P=ParseInput()  # create instance of the parser class to parse the input file
        try:
            P.parse(input_file)
        except Exception, e:
            print "ERROR: Could not parse input file:\n%s"%e
            sys.exit(1)
        try:
            P.check_variables()
        except Exception, e:
            print "ERROR in input variables!\n"
            print e
            sys.exit(1)
        
        try:    #    read in a pdb file and generate appropriate graph (as specified by the input parameters)
            PDB = PDBData(P.pdb_file,P.chains,P.residues,P.distance_type,P.short_distance_threshold,P.long_distance_threshold,P.merge_proRS,P.include_ligand,P.ligand_chain,P.ligand_name,P.ligand_atoms)
        except Exception, e:
            print "Error when creating instance of the class PDBData"
            print e
            sys.exit(1)        
        try: #    read in NMR files and generate appropriate graph (as specified by the input parameters)
            if P.include_ligand:
                NMR = NMRData(P.labelled_residues,P.read_extracted_NOEs,P.extracted_NOE_list,P.protein_ligand_noe_file,P.check_NOE_reciprocity,P.merge_NOEs,P.merge_LV_label,P.include_ligand_noes)
            else:
                NMR = NMRData(P.labelled_residues,P.read_extracted_NOEs,P.extracted_NOE_list,"",P.check_NOE_reciprocity,P.merge_NOEs,P.merge_LV_label)

        except Exception, e:
            print "Error when creating instance of the class NMRData"
            print e
            sys.exit(1)
        ###
        positions,extended_info,info = PDB.extract_positions_indices(P.pdb_flag)
        dist_table_short,info_short,dist_matrix_short,mask_dist_short = PDB.interaction_network(positions,info,P.short_distance_threshold,P.dist_matrix,P.indx_file)
        dist_table_long,info_long,dist_matrix_long,mask_dist_long = PDB.interaction_network(positions,info,P.long_distance_threshold,P.dist_matrix,P.indx_file)
        
        pdb_dict_short = PDB.get_connection_dict(dist_table_short,info_short,False)
        pdb_dict_long = PDB.get_connection_dict(dist_table_long,info_long,False)
        
        if P.include_ligand:
            print "Adding ligand connections to the protein structure graph!"
            lig_dist_table,lig_info,lig_dist_matrix,lig_dist_mask=PDB.ligand_interaction_network(positions,extended_info,P.ligand_distance_threshold)
                
        if P.merge_LV_label:
            print "Leucine and Valine labels are not known\nMerging..."
            pdb_dict_short = PDB.turn_val_to_leu(pdb_dict_short)
            pdb_dict_long = PDB.turn_val_to_leu(pdb_dict_long)
        else:
            print "Leucines and Valines labels are known!"       
        
        if P.include_ligand:
            pdb_dict_lig = PDB.get_connection_dict(lig_dist_table,lig_info,True)
            pdb_dict_short= PDB.join_ligand_connections(pdb_dict_short,pdb_dict_lig)
            pdb_dict_long= PDB.join_ligand_connections(pdb_dict_long,pdb_dict_lig)
            short_node_list,short_adjacency = PDB.get_nodes_adjacency(pdb_dict_short,True,P.merge_proRS)
            long_node_list,long_adjacency = PDB.get_nodes_adjacency(pdb_dict_long,True,P.merge_proRS)
        else:
            short_node_list,short_adjacency = PDB.get_nodes_adjacency(pdb_dict_short,False,P.merge_proRS)
            long_node_list,long_adjacency = PDB.get_nodes_adjacency(pdb_dict_long,False,P.merge_proRS)
        
        #NMR.CheckInputFormats()
        #NMR.InputDimensions()
        NMR.set_noe_dict()
        
        if P.include_ligand_noes:
            noe_dict = NMR.join_ligand_connections(NMR.noes,NMR.lig_noes)
            noe_node_list,noe_adjacency = NMR.get_nodes_adjacency(noe_dict,True,P.merge_proRS)
            return P.labelled_residues,noe_node_list,noe_adjacency,short_node_list,short_adjacency,True
        else:
            noe_node_list,noe_adjacency = NMR.get_nodes_adjacency(NMR.noes,False,P.merge_proRS)
            return P.labelled_residues,noe_node_list,noe_adjacency,short_node_list,short_adjacency,False
        
    def check_graphs(self,data_vertices,data_adjacency,struct_vertices,struct_adjacency,Gp,min_size):
        ##    remove disconnected nodes from the G1 - smaller graph
        data_vertices,data_adjacency = Gp.remove_disconnected_nodes(data_vertices,data_adjacency)
        try:
            # //TO DO//this bugs occasionally! check!!
            data_vertices,data_adjacency = Gp.remove_small_subgraphs(data_vertices,data_adjacency,min_size)
        except Exception,e:
            print "Error in method remove_small_subgraphs (class Graph)%s"%(e)
            print "Proceeding without removal"
            #sys.exit(1)
            
        # basic graphs checks
        small_graph = Gp.networkx_graph(data_vertices,data_adjacency)
        big_graph = Gp.networkx_graph(struct_vertices,struct_adjacency)
        
        Gp.check_graphs_size(small_graph,big_graph)
        Gp.check_graphs_labels(small_graph,big_graph)
        print "Basic NOE graph analysis:\n"
        print "N (nodes, NOE graph) = ",len(small_graph.nodes())
        print "N (edges, NOE graph) = ",len(small_graph.edges())
        print "N (nodes, PDB graph) = ", len(big_graph.nodes())
        print "N (edge, PDB graph) = ", len(big_graph.edges())
        print "NOE sparsity >> ", len(small_graph.edges())/float(len(big_graph.edges()))
        print "average degree >> ",np.mean(nx.degree(small_graph).values())
        print "median degree >>",np.median(nx.degree(small_graph).values())
        
        return data_vertices,data_adjacency,struct_vertices,struct_adjacency
    
    def subgraph_isomorphism(self,noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,Gp,tag):
        
        graph_noe,graph_noe_indexing = Gp.igraph_graph(noe_vertices,noe_adjacency)
        graph_structure,graph_structure_indexing = Gp.igraph_graph(structure_vertices,structure_adjacency)
        EP = IgraphSubIso() #    instance of a class for subgraph isomorphism check and extraction
        start_vf2 = time.time() # measure how long it takes for subgraph isomorphism to be evaluated
        if EP.igraph_subisomorphism(graph_noe,graph_structure):
            n_solutions = EP.igraph_count_subisomorphism(graph_noe,graph_structure)
            print "Found >> ", n_solutions, "SI solutions!"
            all_results = EP.igraph_list_subisomorphism(graph_noe,graph_structure,graph_noe_indexing,graph_structure_indexing,[],[],noe_vertices,[],structure_vertices,[],[],False)
            labels_result_dict=EP.igraph_subisomorphism_write_result(self.outdir+os.sep+"vf2",str(tag),all_results,graph_noe_indexing,graph_structure_indexing)
            #print "SI condition met for the entire NOE graph!\nN (explained NOEs) = ",len(small_graph.edges())
            end_vf2 = time.time()-start_vf2
            print "Took >> ", end_vf2, " sec to evaluate all subgraph isomorphisms"       
            return True
        else:
            return False    # subgraph isomorphism not detected
     
    # @param noe_vertices NMR graph vertices
    # @param noe_adjacency nested list of NMR graph vertices' adjacency relationships
    # @param structure_vertices structure graph vertices
    # @param structure_adjacency nested list of structure graph vertices' adjacency relationships
    # @param Gp an instance of the class Graph (defined in graph_heuristics.py)
    # @param strip_mode a heuristic mode that allows reducing vertex assignment possibilities to predefined 'filter' hash table
    # @filter a hash table - if used it will reduce candidates_dict i.e. vertex assignment options to the values encountered in the filter; by default set to empty 
    # @retval a hash table with nmr graph vertices as keys and structure graph vertices (their assignment options) as values
    def set_assignment_options(self,noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,Gp,strip_mode,tag,filter_dict):
        
        candidates_dict = Gp.compare(noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,"Jaccard",False)
        H = Heuristic() # create an instance of Heuristics
        # execute graph matching order heuristics based on the Jaccard coefficients for direct neighbourhoods and the Hungarian algorithm 
        hungarian_output = self.outdir+os.sep+"hungarian_assignment"+str(tag)+".out"
        munkres_candidates = H.hungarian_first_order(noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,candidates_dict,hungarian_output,1)
        candidates_dict = Gp.combine_munkres_candidates(candidates_dict,munkres_candidates,False) 
        ##
        if strip_mode=='on' and len(filter_dict.keys())>0:   # execute only if the filter is non-empty
            print "Running with strip mode -- reducing assignment options to filter"
            candidates_dict = H.update_dict(candidates_dict,filter)

        return candidates_dict
    
    ## runs short iterations of mcgregor mces algorithm to generate the best starting vertex for the search and the best order of vertices for matching
    # @param Gp instance of the class Graph -- created upon data/parameter parsing
    # @param assign_options hash table of vertex assignment options with nmr-vertices as keys and structure-vertices as values
    # @param noe_vertices NMR graph vertices
    # @param noe_adjacency nested list of NMR graph vertices' adjacency relationships
    # @param structure_vertices structure graph vertices
    # @param structure_adjacency nested list of structure graph vertices' adjacency relationships
    # @retval fin_candidates an optimised order hash table of vertex assignment options with nmr-vertices as keys and structure-vertices as values
    # @retval fin_vertices an optimised order of of data (NMR) graph vertices for the MCES search
    # @retval fin_adjacencies a vertex adjacency order that follows the order established by fin_vertices
    # @retval start_time starting time of the optimisation calculation
    def optimise_run_order(self,Gp,candidates_dict,noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,n_prior_iter,n_prior_mces,max_time,tag,version="py",optimise_mode='y'):               
        
        start_time = time.time() # start a timer
        # create output files for convergence analysis
        #fout1=open(self.outdir+os.sep+'conv.1.',"w")
        #fout2=open(self.outdir+os.sep+'conv.2.',"w")
        #fout1.close()
        #fout2.close()
        # start empty storage data structures
        vertex_lists = {}
        vertex_adjacencies = {}   
        vertex_scores = {}
        final_lists = {}  
        tot_n_iter = 0   
        while tot_n_iter < n_prior_iter:    # run specified number of priority iterations (from the input text file)
            tot_n_iter+=1
            for n in range(len(noe_vertices)):
                # set as the first vertex each of the nmr (data) graph vertices
                shuffled_noe_vertices = Gp.optimal_graph_order_from_node(noe_vertices[n],noe_vertices,noe_adjacency)
                shuffled_noe_adjacency = Gp.re_order_adjecancy(shuffled_noe_vertices,noe_vertices,noe_adjacency)
                """
                //TO CHECK//    deepcopying probably not necessary below -- omit 
                """
                vertex_lists.setdefault(copy.deepcopy(noe_vertices[n]),copy.deepcopy(shuffled_noe_vertices)) 
                vertex_adjacencies.setdefault(copy.deepcopy(noe_vertices[n]),copy.deepcopy(shuffled_noe_adjacency))
                
                if version=="py":
                    MP = MCES_PY(noe_vertices)
                    MP.prepare_mces_py(noe_adjacency,structure_vertices,structure_adjacency,structure_vertices,structure_adjacency,candidates_dict)
                    iter_score,results = MP.mcgregor_main(self.outdir+os.sep+"optimize_mces_step1_py.out",False,1,False,None,False)
                    
                if(iter_score!=0):  # if mces of size>0 is found for this starting vertex -- store the iteration scores and vertices that score them in a dictionary
                    # deepcopying probably not neccessary below                    
                    vertex_scores.setdefault(iter_score,[]).append(copy.deepcopy(noe_vertices[n])) # save all the vertices with this score in the dictionary
                                         
            try:
                best_start_vertex = random.choice(vertex_scores[max(vertex_scores.keys())])    # choose randomly from a pool of the best scoring vertices
                print "Started the mapping search from >> ",best_start_vertex
                sys.stdout.flush()
                
                if version=="py":
                    MP = MCES_PY(vertex_lists[best_start_vertex])
                    MP.prepare_mces_py(vertex_adjacencies[best_start_vertex],structure_vertices,structure_adjacency,structure_vertices,structure_adjacency,candidates_dict)
            
            except ValueError:    #except empty vertex_scores dictionary (i.e. if the first step of the optimization did not generate any mces -- start the search first vertex)
                best_start_vertex = 0
                print "Started the mapping search from >> ",best_start_vertex
                sys.stdout.flush()
                ##         
                shuffled_noe_vertices = Gp.optimal_graph_order_from_node(noe_vertices[0],noe_vertices,noe_adjacency)
                shuffled_noe_adjacency = Gp.re_order_adjecancy(shuffled_noe_vertices,noe_vertices,noe_adjacency)
                
                if version=="py":
                    MP=MCES_PY(shuffled_noe_vertices)
                    MP.prepare_mces_py(shuffled_noe_adjacency,structure_vertices,structure_adjacency,structure_vertices,structure_adjacency,candidates_dict)
                    
            if version=="py":
                iter_score,best_results=MP.mcgregor_main(self.outdir+os.sep+"optimize_mces_step2_py.out",True,n_prior_mces,True,max_time,False)
                            
            #write conv.2.out file -- monitor distribution of vertex scores in the optimization
            fout2=open(self.outdir+os.sep+"conv.2."+tag+".out","a")
            keys=vertex_scores.keys()
            if(len(keys)==0):
                max_score=0
            else:
                max_score=np.max(keys)
            for i in range((max_score+2)):
                if i in keys:
                    fout2.write("%i\t%i\t%i\n" % (tot_n_iter,i,len(set(vertex_scores[i]))))  # write out how many vertices scored this score
                else: # no vertices scored this score
                    fout2.write("%i\t%i\t%i\n" % (tot_n_iter,i,0))
            fout2.write("\n\n")
            fout2.close()

            # monitor scores of mces solutions accummulated during optimization
            fout1=open(self.outdir+os.sep+"conv.1."+tag+".out","a")
            fout1.write("%i\t%i\n"%(tot_n_iter,iter_score)) # write out optimization step (col1) and score of collected mces (col2)
            fout1.close()
            
            print "In iteration ",tot_n_iter,"scored >> ", iter_score
            if (iter_score!=0) and (len(vertex_lists.keys())!=0):
                final_lists[iter_score]=(candidates_dict,vertex_lists[best_start_vertex],vertex_adjacencies[best_start_vertex])
            elif (iter_score!=0) and (len(vertex_lists.keys())==0):
                #// This condition unnecessary// 'full' mode always applies
                final_lists[iter_score]=(candidates_dict,shuffled_noe_vertices,shuffled_noe_adjacency)
            elif iter_score==0:
                print "Iteration led to no MCES "
                #sys.exit(1)

            # after iterating over each starting vertex possibility, reorder candidates_dict according to the best_results   
            
            if version=="py":
                candidates_dict = MP.re_prioritize(best_results,candidates_dict)
            
        #self.write_out_gnu(tag)
        
        if not bool(final_lists): # if final_lists dictionary is empty -- iterations did not lead to MCES -- error somewhere
            print 'No solutions found in any of the optimization runs! An error in settings or input data. Aborting.'
            #return candidates_dict,shuffled_noe_vertices,shuffled_noe_adjacency,start_time
            sys.exit(1)
        
        else:
            best_start_vertex = max(final_lists.keys())    # choose at random start from the pool of the best solutions from all iterations
            fin_candidates=final_lists[best_start_vertex][0]    # final candidates dict
            fin_vertices=final_lists[best_start_vertex][1]  # final vertex order
            fin_adjacencies=final_lists[best_start_vertex][2]   # final vertex adjacency
            return fin_candidates,fin_vertices,fin_adjacencies,start_time
        
        #else:
        #    shuffled_noe_node_list = self.G.optimal_graph_order_from_node(best_start_node,noe_node_list,noe_adjacency)   # final NOE nodes list
        #    shuffled_noe_adjacency = self.G.re_order_adjecancy(shuffled_noe_node_list,noe_node_list,noe_adjacency) #    final adjacency list
        #    return candidates_dict,shuffled_noe_node_list,shuffled_noe_adjacency
        
    ## on the basis of the data structures exiting the optimization step -- executes complete McGregor run
    
    ## run python version of mcgregor mces algorithm
                    
    def run_complete_mcgregor_py(self,candidates_dict,noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,time_elapsed,tag,mces_mode="all"):
        
        MC = MCES_PY(noe_vertices)
        MC.prepare_mces_py(noe_adjacency,structure_vertices,structure_adjacency,structure_vertices,structure_adjacency,candidates_dict)
        
        outfile=self.outdir+os.sep+"mces_"+str(tag)+"_py.txt"
                
        if mces_mode == "one":
            print "\nstarting exact search for a MCES...\n"            
            iter_score,best_results = MC.mcgregor_main(outfile,False,None,False,None,False)
            time_took = time.time() - time_elapsed
            print 'Took >> ', time_took, " to find 1 MCES"
            
        elif mces_mode == "all":
            print "\nstarting exact search for all MCES...\n"
            iter_score,best_results = MC.mcgregor_main(outfile,True,None,False,None,False)
            time_took = time.time() - time_elapsed
            print 'Took >> ', time_took, " to find all MCES"
        else:   # potentially unnecessary - checked by Parser
            print "Unrecognized mcesmode in mcgregor_main!"
            sys.exit(1)
         
                    
    ## prepares and executes a gnuplot script to plot conv1.tag.out and conv.2.tag.out
    # @param runtag - tag that defines conv1.tag.out and conv2.tag.out file names
    def write_out_gnu(self,runtag):
        
        gnu=open(self.outdir+os.sep+'conv.gp','w')
        #gnu.write('set term post eps enh color solid\n')
        gnu.write('set output \'%s\'\n' % (self.outdir+os.sep+'conv'+runtag+'.eps'))
        gnu.write('set size square\n')
        gnu.write('set cblabel \'priority iteration number\'\n')
        gnu.write('set xlabel \'mces size\'\n')
        gnu.write('set ylabel \'count\'\n')
        gnu.write('set key left\n')
        gnu.write('plot \\\n')
        gnu.write('\'%s/conv.1.%s.out\' u 2:1 ti \'progress of best\' w li lw 3,\\\n' % (self.outdir,runtag))
        gnu.write('\'%s/conv.2.%s.out\' u 2:3:1 noti w li lc palette\n' % (self.outdir,runtag))
        gnu.close()
        # assumes gnuplot installed on the system
        os.system('gnuplot '+self.outdir+os.sep+"conv.gp")        
        
    ## splits NOE connections to separate files for separate connected components
    # @ retval noe_subgraphs - list of graph objects of networkx module
    def split_data_graph(self,noe_vertices,noe_adjacency,Gp):
        
        noe_subgraphs = Gp.split_connections_over_subgraphs(noe_vertices,noe_adjacency)
        
        return noe_subgraphs
    
    ## takes graph object of networkx class and returns vertex list and nested adjacency list
    # @ retval subgraph_vertices - list of nmr data subgraph vertices
    # @ retval subgraph_adjacency - nested list of nmr data subgraph vertex connectivity (adjacency)
    def get_subgraph_data(self,conn_subgraph,Gp):
        
        subgraph_conn_dict = Gp.get_conns_dict(conn_subgraph)
        subgraph_vertices,subgraph_adjacency = Gp.get_nodes_adjacency(subgraph_conn_dict)
        
        return subgraph_vertices,subgraph_adjacency

    def analyse_magma_results(self):        
        
        outfile=self.outdir+os.sep+"combined_result.res"
        print 'Writing combined results to: ',outfile
        result_dict=analysis.analyse_all(self.outdir)
        analysis.write_result_dict(result_dict,outfile)
