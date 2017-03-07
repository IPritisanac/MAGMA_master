"""
__author__: Iva Pritisanac

relies on network analysis package igraph | igraph.org
uses igraph implementation of the algorithm VF2 to evaluate, count and write graph-subgraph isomorphism
"""

from generic import GenericMethods
import igraph as ig
#import matplotlib.pyplot as plt
import copy
import numpy as np

##    class prepares the data structures for methods of imported packages
#    calls methods of the imported packages igraph and networkx
class IgraphSubIso(GenericMethods):
            
    def igraph_subisomorphism(self,small_graph,big_graph):        
        return big_graph.subisomorphic_vf2(small_graph,color1 =big_graph.vs["color"] ,color2 =small_graph.vs["color"])
    
    def igraph_count_subisomorphism(self,small_graph,big_graph):    
        N_subisomorphisms = big_graph.count_subisomorphisms_vf2(small_graph,color1 = big_graph.vs["color"] ,color2 = small_graph.vs["color"])
        return N_subisomorphisms
     
    
    def create_result_dict(self,subisomorphisms_list):
        #    index in the sublist is index1
        results = {}
        for list in subisomorphisms_list:
            for e in range(len(list)):
                results.setdefault(e,[]).append(list[e])
        for key,value in results.items():
            results[key] = set(value)
        return results
        
    def igraph_list_subisomorphism(self,small_graph,big_graph,index1,index2,edge_matrix1,edge_matrix2,noe_nodes,noes,pdb_nodes,pdbs,pdbl,rescore):
        
        total_results_score = {}
        all_isomorphisms = big_graph.get_subisomorphisms_vf2(small_graph,color1=big_graph.vs["color"],color2=small_graph.vs["color"])
        
        if rescore == True:     
            for list in all_isomorphisms:
                assignment = {index1[l]:index2[list[l]] for l in range(len(list))}
                mat1 = copy.deepcopy(edge_matrix1)
                mat2 = copy.deepcopy(edge_matrix2)
                score_short = self.SetEdgeAssignment(mat1,assignment,noe_nodes,noes,pdb_nodes,pdbs)
                score_long = self.SetEdgeAssignment(mat2,assignment,noe_nodes,noes,pdb_nodes,pdbl)
                total_results_score = self.AssembleResults(total_results_score,score_short,score_long,assignment)
            return total_results_score
        else:
            results = self.create_result_dict(all_isomorphisms)
            return results
        return total_results_score
    
    def igraph_subisomorphism_write_result(self,output_file,index_file,assign_dict,init_index_small,init_index_big):
        label_results = {} 
        out = open(output_file+"_"+index_file+".txt","w")        
        for methyl,assignments in assign_dict.items():
            for assign in assignments:
                out.write("%s\t%s\n"%(init_index_small[methyl],init_index_big[assign]))
                label_results.setdefault(init_index_small[methyl],[]).append(init_index_big[assign])
        out.close()
        return label_results
    
    