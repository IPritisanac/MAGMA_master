#! /usr/bin/python

"""
executes methods defined in mces.py (c) Iva Pritisanac | iva.pritisanac@gmail.com
"""

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

class Magma():
    
    def __init__(self,inputf):
        
        self.dirpath=os.path.dirname(os.path.abspath(inputf)) # extract directory name out of the input filename
        if not os.path.exists(self.dirpath+os.sep+"results"): # if the output results directory does not exist on this path
            os.makedirs(self.dirpath+os.sep+"results")   # create the directory
        self.outdir=self.dirpath+os.sep+"results"   # set the output directory
        
    ## run parser and return only variables -- different modes, integer parameter values etc.    
    def parse_variables(self,input_file):
        
        P=ParseInput()  # create instance of the parser class to parse the input file
        try:
            P.parse(input_file)
        except Exception, e:
            print "Could not parse input file:\n%s"%e
            sys.exit(1)
        try:
            P.check_variables()
        except Exception, e:
            print "ERROR in input variables!\n"
            print e
            sys.exit(1)
            
        return P.mode,P.min_size,P.mces_mode,P.strip_mode,P.prior_iter,P.n_prior_mces,P.max_run_time
    
    ## run parser to fill out and return data structures
    def parse_data(self,input_file):

        P=ParseInput()  # create instance of the parser class to parse the input file
        try:
            P.parse(input_file)
        except Exception, e:
            print "Could not parse input file:\n%s"%e
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
                NMR = NMRData(P.HMQC_list,P.format_peak_list,P.labelled_residues,P.read_NOEs_from_file,P.NOE_list,P.format_NOE_list,P.dimensions_NOEs,P.read_extracted_NOEs,P.extracted_NOE_list,P.protein_ligand_noe_file,P.check_NOE_reciprocity,P.check_NOE_intensity,P.merge_NOEs,P.merge_LV_label,P.include_ligand_noes)
            else:
                NMR = NMRData(P.HMQC_list,P.format_peak_list,P.labelled_residues,P.read_NOEs_from_file,P.NOE_list,P.format_NOE_list,P.dimensions_NOEs,P.read_extracted_NOEs,P.extracted_NOE_list,"",P.check_NOE_reciprocity,P.check_NOE_intensity,P.merge_NOEs,P.merge_LV_label)

        except Exception, e:
            print "Error when creating instance of the class NMRData"
            print e
            sys.exit(1)
        ###
        positions,extended_info,info = PDB.ExtractPositionsIndices(P.pdb_flag)
        dist_table_short,info_short,dist_matrix_short,mask_dist_short = PDB.InteractionNetwork(positions,info,P.short_distance_threshold,P.dist_matrix,P.indx_file)
        dist_table_long,info_long,dist_matrix_long,mask_dist_long = PDB.InteractionNetwork(positions,info,P.long_distance_threshold,P.dist_matrix,P.indx_file)
        
        pdb_dict_short = PDB.GetConnectionDict(dist_table_short,info_short,False)
        pdb_dict_long = PDB.GetConnectionDict(dist_table_long,info_long,False)
        
        if P.include_ligand:
            print "Adding ligand connections to the protein structure graph!"
            lig_dist_table,lig_info,lig_dist_matrix,lig_dist_mask=PDB.LigandInteractionNetwork(positions,extended_info,P.ligand_distance_threshold)
        
        """
        //TO DO// move this to the plotting scripts
        
        #PDB.WriteDistanceFile(pdb_dict_long_distances,"dist.txt")
        #PDB.PrettyPlotContacts(dist_table_short, positions, dist_matrix_short, mask_dist_short, P.short_distance_threshold)
        #PDB.PrettyPlotContacts(dist_table_long, positions, dist_matrix_long, mask_dist_long, P.long_distance_threshold)
        ## if we do not discriminate Leucines from Valines and vice versa
        """
        
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
            short_node_list,short_adjacency = PDB.GetNodesAdjacency(pdb_dict_short,True,P.merge_proRS)
            long_node_list,long_adjacency = PDB.GetNodesAdjacency(pdb_dict_long,True,P.merge_proRS)
        else:
            short_node_list,short_adjacency = PDB.GetNodesAdjacency(pdb_dict_short,False,P.merge_proRS)
            long_node_list,long_adjacency = PDB.GetNodesAdjacency(pdb_dict_long,False,P.merge_proRS)
        
        NMR.CheckInputFormats()
        NMR.InputDimensions()
        NMR.SetNOEDict()
        
        if P.include_ligand_noes:
            noe_dict = NMR.JoinLigandConnections(NMR.noes,NMR.lig_noes)
            noe_node_list,noe_adjacency = NMR.GetNodesAdjacency(noe_dict,True,P.merge_proRS)
            return P.labelled_residues,noe_node_list,noe_adjacency,short_node_list,short_adjacency,True
        else:
            noe_node_list,noe_adjacency = NMR.GetNodesAdjacency(NMR.noes,False,P.merge_proRS)
            return P.labelled_residues,noe_node_list,noe_adjacency,short_node_list,short_adjacency,False
        
    def check_graphs(self,data_vertices,data_adjacency,struct_vertices,struct_adjacency,Gp,min_size):
        
        ##    remove disconnected nodes from the G1 - smaller graph
        noe_node_list,noe_adjacency = Gp.RemoveDisconnectedNodes(data_vertices,data_adjacency)
        try:
            """
            This bugs on ubiqutin! Check!!
            """
            data_vertices,data_adjacency = Gp.RemoveSmallSubgraphs(data_vertices,data_adjacency,min_size)
        except Exception,e:
            print e
            sys.exit(1)
            
        # basic graphs checks
        small_graph = Gp.NetworkxGraph(data_vertices,data_adjacency)
        big_graph = Gp.NetworkxGraph(struct_vertices,struct_adjacency)
        
        """
        //TO DO// move this to the plotting scripts
        #G.PlotGraphLabels(small_graph,"example.pdf")
        """
        Gp.CheckGraphsSize(small_graph,big_graph)
        Gp.CheckGraphsLabels(small_graph,big_graph)
        print "Basic NOE graph analysis:\n"
        print "N (nodes, NOE graph) = ",len(small_graph.nodes())
        print "N (edges, NOE graph) = ",len(small_graph.edges())
        print "N (nodes, PDB graph) = ", len(big_graph.nodes())
        print "N (edge, PDB graph) = ", len(big_graph.edges())
        print "NOE sparsity >> ", len(small_graph.edges())/float(len(big_graph.edges()))
        print "average degree >> ",np.mean(nx.degree(small_graph).values())
        print "median degree >>",np.median(nx.degree(small_graph).values())
        
        return data_vertices,data_adjacency,struct_vertices,struct_adjacency
    
    def subgraph_isomorphism(self,noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,tag,Gp):
        
        graph_noe,graph_noe_indexing = Gp.IgraphGraph(noe_vertices,noe_adjacency)
        graph_structure,graph_structure_indexing = Gp.IgraphGraph(structure_vertices,structure_adjacency)
        EP = IgraphSubIso() #    instance of a class for subgraph isomorphism check and extraction
        start_vf2 = time.time() # measure how long it takes for subgraph isomorphism to be evaluated
        if EP.IgraphSubIsomorphism(graph_noe,graph_structure):
            n_solutions = EP.IgraphCountSubIsomorphism(graph_noe,graph_structure)
            print n_solutions
            all_results = EP.IgraphListSubIsomorphism(graph_noe,graph_structure,graph_noe_indexing,graph_structure_indexing,[],[],noe_vertices,[],structure_vertices,[],[],False)
            labels_result_dict=EP.IgraphSubIsomorphismWriteResult(self.outdir+os.sep+"vf2",str(tag),all_results,graph_noe_indexing,graph_structure_indexing)
            #print "SI condition met for the entire NOE graph!\nN (explained NOEs) = ",len(small_graph.edges())
            end_vf2 = time.time()-start_vf2
            print "Took >> ", end_vf2, " sec to evaluate all subgraph isomorphisms"       
            return True
        else:
            return False    # subgraph isomorphism not detected
        
    ## run magma mcgregor algorithm
    #@param noe_vertices NMR graph vertices
    #@param noe_adjacency nested list of NMR graph vertices' adjacency relationships
    #@param structure_vertices structure graph vertices
    #@param noe_adjacency nested list of structure graph vertices' adjacency relationships
    #@param Gp an instance of the class Graph (defined in graph_heuristics.py)
    #@param mces_mode a mode that defines if mcgregor algorithm is to look for 'one' or 'all' mces solutions
    #@param strip_mode a heuristic mode that allows reducing vertex assignment possibilities to predefined 'filter' hash table
    #@param order_mode a mode that defines optimization of running order ??
    #@param tag a tag added to ??
    #@filter a hash table - if used it will reduce candidates_dict i.e. vertex assignment options to the values encountered in the filter; by default set to empty
    #@retval
    #
    def magma_mcgregor(self): #,noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,Gp,mces_mode,strip_mode,order_mode="full",tag='',filter_dict={}):
        
        time_start=time.time() # start measuring run time
        #self.optimise_run_order(candidates_dict,order_mode,noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,20,3,100,tag,'y')
     
    # @param noe_vertices NMR graph vertices
    # @param noe_adjacency nested list of NMR graph vertices' adjacency relationships
    # @param structure_vertices structure graph vertices
    # @param structure_adjacency nested list of structure graph vertices' adjacency relationships
    # @param Gp an instance of the class Graph (defined in graph_heuristics.py)
    # @param strip_mode a heuristic mode that allows reducing vertex assignment possibilities to predefined 'filter' hash table
    # @filter a hash table - if used it will reduce candidates_dict i.e. vertex assignment options to the values encountered in the filter; by default set to empty 
    # @retval a hash table with nmr graph vertices as keys and structure graph vertices (their assignment options) as values
    def set_assignment_options(self,noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,Gp,strip_mode,filter_dict):
        
        candidates_dict = Gp.Compare(noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,"Jaccard",False)
        H = Heuristic() # create an instance of Heuristics
        # execute graph matching order heuristics based on the Jaccard coefficients for direct neighbourhoods and the Hungarian algorithm 
        munkres_candidates = H.HungarianFirstOrder(noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,candidates_dict,1)
        candidates_dict = Gp.CombineMunkresCandidates(candidates_dict,munkres_candidates,False) 
        ##
        if strip_mode=='on' and len(filter_dict.keys())>0:   # execute only if the filter is non-empty
            print "Running with strip mode -- reducing assignment options to filter"
            candidates_dict = H.CraicD(candidates_dict,filter)
            
        return candidates_dict
    
    ## runs short iterations of mcgregor mces algorithm to generate the best starting vertex for the search and the best order of vertices for matching
    # @param Gp instance of the class Graph -- created upon data/parameter parsing
    # @param assign_options hash table of vertex assignment options with nmr-vertices as keys and structure-vertices as values
    # @param noe_vertices NMR graph vertices
    # @param noe_adjacency nested list of NMR graph vertices' adjacency relationships
    # @param structure_vertices structure graph vertices
    # @param structure_adjacency nested list of structure graph vertices' adjacency relationships
    # @retval ??        
    def optimise_run_order(self,Gp,candidates_dict,noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,n_prior_iter,n_prior_mces,max_time,tag,version="c",optimise_mode='y'):               
        
        start_time = time.time() # start a timer
        # create output files for convergence analysis
        fout1=open(self.outdir+os.sep+'conv.1',"w")
        fout2=open(self.outdir+os.sep+'conv.2',"w")
        fout1.close()
        fout2.close()
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
                shuffled_noe_vertices = Gp.OptimalGraphOrderFromNode(noe_vertices[n],noe_vertices,noe_adjacency)
                shuffled_noe_adjacency = Gp.ReOrderAdjecancy(shuffled_noe_vertices,noe_vertices,noe_adjacency)
                """
                //TO CHECK//    deepcopying probably not neccessary below -- omit 
                """
                vertex_lists.setdefault(copy.deepcopy(noe_vertices[n]),copy.deepcopy(shuffled_noe_vertices)) 
                vertex_adjacencies.setdefault(copy.deepcopy(noe_vertices[n]),copy.deepcopy(shuffled_noe_adjacency))
                if version=="c": # run c version of McGregor algorithm
                    MC=McGregor(shuffled_noe_vertices) # instantiate McGregor class given the current order of nmr graph vertices (shuffled_noe_vertices)
                    MC.PrepareMcGregor(shuffled_noe_adjacency,structure_vertices,structure_adjacency,structure_vertices,structure_adjacency,candidates_dict,self.outdir)
                    if optimise_mode=="y":  # 04/03/2017 -- the earlier name 'vanilla' mode
                        iter_score,results = MC.McGregorLoop(False,False,1,False,max_time,False,0)
                    else:
                        iter_score,results = MC.McGregorLoop(False,True,n_prior_mces,True,max_time,False,0)
                if version=="py":
                    pass
                    """
                    // TO DO//
                    """
                if(iter_score!=0):  # if mces of size>0 is found for this starting vertex -- store the iteration scores and vertices that score them in a dictionary
                    """
                    //TO CHECK// deepcopying probably not neccessary here
                    """
                    vertex_scores.setdefault(iter_score,[]).append(copy.deepcopy(noe_vertices[n])) # save all the vertices with this score in the dictionary
                                         
            try:
                best_start_vertex = random.choice(vertex_scores[max(vertex_scores.keys())])    # choose randomly from a pool of the best scoring vertices
                print "Started the mapping search from >> ",best_start_vertex
                sys.stdout.flush()
                #makeInputfile for mcesCore
                MC=McGregor(vertex_lists[best_start_vertex])    # start the search from this vertex, the bfs from that vertex, and adjacency ordered accordingly
                MC.PrepareMcGregor(vertex_adjacencies[best_start_vertex],structure_vertices,structure_adjacency,structure_vertices,structure_adjacency,candidates_dict,self.outdir)
            except ValueError:    #except empty vertex_scores dictionary (i.e. if the first step of the optimization did not generate any mces -- start the search first vertex)
                best_start_vertex = 0
                print "Started the mapping search from >> ",best_start_vertex
                sys.stdout.flush()
                ##         
                shuffled_noe_vertices = Gp.OptimalGraphOrderFromNode(noe_vertices[0],noe_vertices,noe_adjacency)
                shuffled_noe_adjacency = Gp.ReOrderAdjecancy(shuffled_noe_vertices,noe_vertices,noe_adjacency)
                MC=McGregor(shuffled_noe_vertices)
                MC.PrepareMcGregor(shuffled_noe_adjacency,structure_vertices,structure_adjacency,structure_vertices,structure_adjacency,candidates_dict,self.outdir)
                
            if(optimise_mode=="y"):
                iter_score,best_results = MC.McGregorLoop(False,True,n_prior_mces,True,max_time,False,0)
            else:
                iter_score,best_results = MC.McGregorLoop(False,True,n_prior_mces,True,10*max_time,False,0)
            
            #write conv.2.out file -- monitor distribution of vertex scores in the optimization
            fout2=open(self.outdir+os.sep+"conv.2"+tag+".out","a")
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
            fout1=open(self.outdir+os.sep+"conv.1"+tag+".out","a")
            fout1.write("%i\t%i\n"%(tot_n_iter,iter_score)) # write out optimization step (col1) and score of collected mces (col2)
            fout1.close()
            
            print "In iteration ",tot_n_iter,"scored >> ", iter_score
            if (iter_score!=0) and (len(vertex_lists.keys())!=0):
                final_lists[iter_score]=(candidates_dict,vertex_lists[best_start_vertex],vertex_adjacencies[best_start_vertex])
            elif (iter_score!=0) and (len(vertex_lists.keys())==0):
                """
                // This condition is unnecessary ?! // 'full' mode always applies
                """
                final_lists[iter_score]=(candidates_dict,shuffled_noe_vertices,shuffled_noe_adjacency)
            elif iter_score==0:
                print "Iteration led to no MCES "
                #sys.exit(1)

            # after iterating over each starting vertex possibility, reorder candidates_dict according to the best_results   
            candidates_dict = MC.RePrioritize(best_results,candidates_dict)
                
        self.write_out_gnu(tag)
        
        if not bool(final_lists): # if final_lists dictionary is empty -- iterations did not lead to MCES -- error somewhere
            print 'No solutions found in any of the optimization runs! An error in settings or input data. Aborting.'
            """
            ABORT HERE!! Do not return anything
            """
            return candidates_dict,shuffled_noe_vertices,shuffled_noe_adjacency,start_time
            #sys.exit(1)
        
        else:
            best_start_vertex = max(final_lists.keys())    # choose at random start from the pool of the best solutions from all iterations
            fin_candidates=final_lists[best_start_vertex][0]    # final candidates dict
            fin_vertices=final_lists[best_start_vertex][1]  # final vertex order
            fin_adjacencies=final_lists[best_start_vertex][2]   # final vertex adjacency
            return fin_candidates,fin_vertices,fin_adjacencies,start_time
        
        #else:
        #    shuffled_noe_node_list = self.G.OptimalGraphOrderFromNode(best_start_node,noe_node_list,noe_adjacency)   # final NOE nodes list
        #    shuffled_noe_adjacency = self.G.ReOrderAdjecancy(shuffled_noe_node_list,noe_node_list,noe_adjacency) #    final adjacency list
        #    return candidates_dict,shuffled_noe_node_list,shuffled_noe_adjacency
        
    ## on the basis of the data structures exiting the optimization step -- executes complete McGregor run
            
    def run_complete_mcgregor(self,candidates_dict,noe_vertices,noe_adjacency,structure_vertices,structure_adjacency,time_elapsed,tag,mces_mode="all"):
        
        MC = McGregor(noe_vertices)    # create an instance of McGregor class with the current list of vertices
        MC.PrepareMcGregor(noe_adjacency,structure_vertices,structure_adjacency,structure_vertices,structure_adjacency,candidates_dict,self.outdir)
        # define name of the output file
        outfile=self.outdir+os.sep+"mces_"+str(tag)+".txt"
        
        if mces_mode=="one":
            iter_score,best_results=MC.McGregorLoop(outfile,False,None,False,None,False,1)
            total_time = time.time()-time_elapsed
            print "Calculation took >> ", total_time, " to find 1 mces"
        elif mces_mode=="all":
            iter_score,best_results=MC.McGregorLoop(outfile,True,None,False,None,False,1)
            total_time = time.time()-time_elapsed
            print "Calculation took >> ", total_time, " to find all mces"
        else:   # potentially unnecessary - check for by Parser
            print "Unrecognized mcesmode in run_complete_mcgregor!"
            sys.exit(1)
                    
    ## prepares and executes a gnuplot script to plot conv1.tag.out and conv.2.tag.out
    # @param runtag - tag that defines conv1.tag.out and conv2.tag.out file names
    def write_out_gnu(self,runtag):
        
        gnu=open(self.outdir+os.sep+'conv.gp','w')
        gnu.write('set term post eps enh color solid\n')
        gnu.write('set output \'%s\'\n' % (self.outdir+'/conv'+runtag+'.eps'))
        gnu.write('set size square\n')
        gnu.write('set cblabel \'priority iteration number\'\n')
        gnu.write('set xlabel \'mces size\'\n')
        gnu.write('set ylabel \'count\'\n')
        gnu.write('set key left\n')
        gnu.write('plot \\\n')
        gnu.write('\'%s/conv.1%s.out\' u 2:1 ti \'progress of best\' w li lw 3,\\\n' % (self.outdir,runtag))
        gnu.write('\'%s/conv.2%s.out\' u 2:3:1 noti w li lc palette\n' % (self.outdir,runtag))
        gnu.close()
        # assumes gnuplot installed on the system
        os.system('gnuplot '+self.outdir+os.sep+"conv.gp")        
        
# Splits NOE connections to separate files for separate connected components

"""
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
    
    #Call implementation of Vf2 algorithm (c++) implemented in python package igraph
    
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
        
        #Apply heuristic before MCES search
    
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
"""
//To do//
Analyse solutions
Plot solutions
"""