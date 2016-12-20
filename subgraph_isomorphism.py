from generic import GenericMethods
import igraph as ig
import matplotlib.pyplot as plt
import copy
import numpy as np

class IgraphSubIso(GenericMethods):
    #    class prepares the data structures for methods of imported packages
    #    calls methods of the imported packages igraph and networkx
    def QuickPlotIsomorphism(self,results_dict,dim1,dim2):
        matrix=np.zeros((dim2,dim2),dtype=np.int)
        xlabels = []
        for key,value in results_dict.items():
            xlabels.append(key)
            for val in value:#VALUE is ATOM
                matrix[key,val] = 1    # rows are peak identifiers

        plt.pcolor(matrix,cmap='Greys')
        plt.xticks(range(0,len(results_dict.keys())+1), fontsize = 10)
        plt.yticks(range(0,len(results_dict.keys())+1), fontsize = 10)
        
        plt.colorbar()
        plt.title('Peak_residue mapping')
        plt.xlabel('Peak ID')
        plt.ylabel('Atom ID')
        plt.show()     
        
    def IgraphSubIsomorphism(self,small_graph,big_graph):        
        return big_graph.subisomorphic_vf2(small_graph,color1 =big_graph.vs["color"] ,color2 =small_graph.vs["color"])
    
    def IgraphCountSubIsomorphism(self,small_graph,big_graph):    
        N_subisomorphisms = big_graph.count_subisomorphisms_vf2(small_graph,color1 = big_graph.vs["color"] ,color2 = small_graph.vs["color"])
        return N_subisomorphisms
    
    def GetIntEdges(self,G1node,G2node,G1nodes,G2nodes,G1edges,G2edges):
        #    information from initial edges-matrix could be used, without the need to check types!
        #    emptied every time the function is called
        #    get all edges of G1 node
        edges1 = []     
        edges2 = []
        indG1=G1nodes.index(G1node)
        indG2=G2nodes.index(G2node)
        for key,value in G1edges.items():
            for val in value[0]:
                if val == indG1:
                    edges1.append(key)
        for key,value in G2edges.items():
            for val in value[0]:
                if val == indG2:
                    edges2.append(key) 
        return edges1,edges2 
    
    def CreateResultDict(self,subisomorphisms_list):
        #    index in the sublist is index1
        results = {}
        for list in subisomorphisms_list:
            for e in range(len(list)):
                results.setdefault(e,[]).append(list[e])
        for key,value in results.items():
            results[key] = set(value)
        return results
        
    def IgraphListSubIsomorphism(self,small_graph,big_graph,index1,index2,edge_matrix1,edge_matrix2,noe_nodes,noes,pdb_nodes,pdbs,pdbl,rescore):
        
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
            results = self.CreateResultDict(all_isomorphisms)
            return results
        return total_results_score
    
    def IgraphSubIsomorphismWriteResult(self,output_file,index_file,assign_dict,init_index_small,init_index_big):
        label_results = {} 
        out = open(output_file+"_"+index_file+".txt","w")        
        for methyl,assignments in assign_dict.items():
            for assign in assignments:
                out.write("%s\t%s\n"%(init_index_small[methyl],init_index_big[assign]))
                label_results.setdefault(init_index_small[methyl],[]).append(init_index_big[assign])
        out.close()
        return label_results
    
    