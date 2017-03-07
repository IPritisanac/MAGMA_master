"""
__author__: Iva Pritisanac 

definition and manipulation of graphs (protein structure graphs and NMR data graphs);
heuristics for ordering of the graph matching search;
munkres algorithm, variation of the Hungarian algorithm (call to the implementation given by the python package: munkres1.0.8 | https://pypi.python.org/pypi/munkres/);
algorithm for subgraph isomorphism (call to the implementation given by the high performance graph library igraph http://igraph.org/python/);
"""

from generic import GenericMethods
import networkx as nx
import random
import sys
import numpy as np
import munkres
import igraph as ig
from munkres import Munkres, print_matrix
import matplotlib.pyplot as plt
from subgraph_isomorphism import IgraphSubIso

class Graph(GenericMethods): 
    #    class to prepare the graphs according to the input parameters
    #    several graph objects can be created    -- for instance, graph objects of packages networkx as igraph
    #    the methods of these respective packages are supported
    
    def __init__(self,residue_types,labels,label_flag="type",ligand=False):
        
        self.label_flag = label_flag
        
        if self.label_flag == "type":
            self.igraph_color_codes = {residue_types[i][0]:i for i in range(len(residue_types))}            
        
        if self.label_flag == "assigned":
            self.igraph_color_codes = {}
            new_label = max(labels.values())+1  # set new label for residue type
            residue_type_color_codes = {residue_types[i][0]:new_label+i for i in range(len(residue_types))}
            for key,value in labels.items():
                self.igraph_color_codes.setdefault(key,value)
            for key,value in residue_type_color_codes.items():
                self.igraph_color_codes.setdefault(key,value)   
                    
        if self.label_flag == "type" and ligand:
            # if ligand atoms are added to the graph
            # add an additional color label
            self.igraph_color_codes["X"]=len(residue_types)+1
            
    ##    method to remove isolated nodes from NOE graph
    #    isolated nodes are those nodes (methyls) in the NOE graph which are not connected through an edge (NOE) to any other nodes (methyls)
    #    removal done both in the node list and in the corresponding adjacency simultaneously                        
    def remove_disconnected_nodes(self,nodes,adjecancy):
       
        nodes_adj_dict = {nodes[n1]:set(adjecancy[n1]) for n1 in range(len(nodes))}        
        node_orphans = [nodes[x] for x in range(len(nodes)) if len(adjecancy[x]) == 0]
        new_node_list = [nodes[n] for n in range(len(nodes)) if nodes[n] not in node_orphans]
        new_adjecancy = [list(nodes_adj_dict[node]) for node in new_node_list]
        
        return new_node_list,new_adjecancy
    
    def reverse_dictionary(self,dictionary):
        if len(dictionary.keys()) == len(dictionary.values()):
            new_dictionary = {value:key for key,value in dictionary.items()}
            return new_dictionary
        else:
            raise Exception("Dictionary to reverse has different number of keys and values!")
    
    def remove_small_subgraphs(self,nodes,adjacency,minimal_size):
        # creates graph object of the networkx 
        # lists connected subgraphs
        # removed those that are smaller than the required minimal size
        nxgraph = self.networkx_graph(nodes,adjacency)
        connected_subgraphs = nx.connected_component_subgraphs(nxgraph)
        subgraphs_dictionary = self.classify_components(connected_subgraphs)       
        updated_nodes,updated_adjacency = self.remove_components(subgraphs_dictionary,nodes,adjacency,minimal_size)        
        return updated_nodes,updated_adjacency      
            
    ##    given the minimum size (min_size) in number of nodes extracts only components > that the minimum     
    def remove_components(self,subgraphs,node_list,adjacency_list,minimum):
        
        remove_nodes = []
        for size,subs in subgraphs.items():
            #print size,subs
            if size <= minimum:
                for subgraph in subs:
                    remove_nodes.extend(subgraph.nodes())
        
        print "Removing small graph components (N(nodes)) <= ", minimum
        print "..."
        nodes_adj_dict = {node_list[n1]:set(adjacency_list[n1]) for n1 in range(len(node_list))}
        new_node_list = [node for node in node_list if node not in remove_nodes]
        new_adjacency = [list(nodes_adj_dict[node]) for node in new_node_list] 
        
        return new_node_list,new_adjacency                   
    
    
    def create_real_val_assign_matrix(self,candidates_dict):
        N_keys = len(candidates_dict.keys())
        #sorted_keys = sorted(candidates_dict.keys(),key = lambda x:int(x[:-1]))
        sorted_keys = sorted(candidates_dict.keys()) #,key = lambda x:int(x[:-1]))
        unique_values = set([v[0] for val in candidates_dict.values() for v in val])
        N_values = len(unique_values)
        #sorted_values = sorted(unique_values,key = lambda x:int(x[:-1]))
        sorted_values = sorted(unique_values) #,key = lambda x:int(x[:-1]))
        assign_matrix = np.zeros(shape=(N_keys,N_values),dtype=int)
        keys_indices = {sorted_keys[k]:k for k in range(len(sorted_keys))}
        values_indices = {sorted_values[k]:k for k in range(len(sorted_values))}
        for key,value in candidates_dict.items():
            for v in value:
                assign_matrix[keys_indices[key],values_indices[v[0]]] = v[1]
        return assign_matrix,keys_indices,values_indices
                
    def munkres_maximise_assignment(self,matrix):
    # call to Munkres algorithm (a variant of the Hungarian algorithm)
    # maximise "profit" along rows of a matrix
        cost_matrix = munkres.make_cost_matrix(matrix,lambda cost: sys.maxint - cost)
        m = Munkres()
        indexes = m.compute(cost_matrix)
        #print_matrix(matrix, msg='Highest profit through this matrix:')
        useful_indices = []
        total = 0
        for row, column in indexes:
            value = matrix[row][column]
            total += value
            #print '(%d, %d) -> %d' % (row, column, value)
            useful_indices.append(((row, column)))
        print 'total score (hungarian) =%d' % total
        return total,useful_indices    
    
    def combine_munkres_candidates(self,dictionary_total,dictionary_munkres,quantitative = False):
        #    for each atom list
        #    create dictionary of values
        #    sort in ascending order
        #    remove munkres/hungarian candidate
        #    give score in range (0, len of ascending order)   
        #    add munkres/hungarian candidate with score *score max +1
        if quantitative:           
            quantitative_candidates_dict = {}
            for peak,atoms in dictionary_total.items():    
                sorted_atoms = {}
                scored_atoms = {}
                for (atom,score) in atoms:
                    if atom!=dictionary_munkres[peak]:
                        sorted_atoms.setdefault(score,[]).append(atom)  # avoid adding candidate from hungarian assignment twice!
                    else:
                        continue                
                scores = sorted(sorted_atoms.keys())
                for s in range(len(scores)):
                    scored_atoms.setdefault(s,sorted_atoms[scores[s]])  # give integer values of scores for eack peak:atom assignment possibility
                max_score = max(scored_atoms.keys())
                scored_atoms[max_score+1] = [dictionary_munkres[peak]]  # give highest value of score to hungarian assignment
                quantitative_candidates_dict.setdefault(peak,scored_atoms)                   
            #print quantitative_candidates_dict
            return quantitative_candidates_dict
           
        else:
            for key,value in dictionary_total.items():
                if key in dictionary_munkres.keys():
                    #print value
                    if dictionary_munkres[key] in value:
                        indx = value.index(dictionary_munkres[key])
                        value.pop(indx) #remove the duplicate
                    value.insert(0,dictionary_munkres[key])
                    dictionary_total[key] = value    #remove duplicates               

            return dictionary_total
    
    def munkres_based_candidates(self,indices,candidates_matrix,indices_rows,indices_columns):
        
        
        munkres_candidates = {}
        reverse_rows = self.reverse_dictionary(indices_rows)
        reverse_columns = self.reverse_dictionary(indices_columns)
        
        for index in indices:
            indx_row = index[0]
            indx_column = index[1]
            print reverse_rows[indx_row],">>",reverse_columns[indx_column]
            munkres_candidates.setdefault(reverse_rows[indx_row],reverse_columns[indx_column])
         
        return munkres_candidates    
        
    def classify_components(self,subgraphs):
        
        subgraphs_dict = {}            
        #for s in range(len(subgraphs)):
        for subgraph in subgraphs:
            subgraphs_dict.setdefault(len(subgraph.nodes()),[]).append(subgraph)
            #subgraphs_dict.setdefault(len(subgraphs[s].nodes()),[]).append(subgraphs[s])
        return subgraphs_dict   
                                        
    def networkx_graph(self,node_list,adjacency_list):        
        graph = nx.Graph()
        graph_edges = [(node_list[i],e) for i in range(len(node_list)) for e in adjacency_list[i]]
        # if unconnected nodes need to be kept in
        # this should be explicitly coded    --> line below
        graph.add_nodes_from(node_list)
        graph.add_edges_from(graph_edges)   # will exclude any unconnected node (orphan) from graph nodes
        return graph
    
    
    def nx_get_components_as_subgraphs(self,graph):
        connected_subgraphs = nx.connected_component_subgraphs(graph)
        return connected_subgraphs

    def split_connections_over_subgraphs(self,nodes,adjacency,file_names = "splitted_connected_graph_"): 
        nx_graph = self.networkx_graph(nodes,adjacency)
        conn_subgraphs = self.nx_get_components_as_subgraphs(nx_graph)        
        
        #for g in range(len(conn_subgraphs)):
        #    subgraph_file = open(file_names + str(g)+".txt","w")
        #    for edge in conn_subgraphs[g].edges():
        #        subgraph_file.write("%s\t%s\n"%(edge[0][-1]+edge[0][:-1]+"C-H",edge[1][-1]+edge[1][:-1]+"C-H"))
        #        subgraph_file.write("%s\t%s\n"%(edge[1][-1]+edge[1][:-1]+"C-H",edge[0][-1]+edge[0][:-1]+"C-H"))
        #    subgraph_file.close()

        return conn_subgraphs
    
    def plot_graph_labels(self,G,figname):
        import matplotlib.pyplot as plt
        plt.figure()
        pos=nx.spring_layout(G) # positions for all nodes
        pos= nx.random_layout(G)
        node_colors = {'I':'y','L':'r','V':'b','A':'m','M':'c','T':'w','X':'k'}
        
        ile_nodes = [node for node in G.nodes() if node[-1] == 'I']
        leu_nodes = [node for node in G.nodes() if node[-1] == 'L']
        val_nodes = [node for node in G.nodes() if node[-1] == 'V']
        ala_nodes = [node for node in G.nodes() if node[-1] == 'A']
        met_nodes = [node for node in G.nodes() if node[-1] == 'M']
        thr_nodes = [node for node in G.nodes() if node[-1] == 'T']
        lig_nodes = [node for node in G.nodes() if node[-1] == 'X']
        
        # nodes
        nx.draw_networkx_nodes(G,pos,
                               nodelist=ile_nodes,
                               node_color=node_colors['I'],
                               node_size=500,
                               alpha=0.8)
        nx.draw_networkx_nodes(G,pos,
                               nodelist=val_nodes,
                               node_color=node_colors['V'],
                               node_size=500,
                               alpha=0.8)
        nx.draw_networkx_nodes(G,pos,
                               nodelist=leu_nodes,
                               node_color=node_colors['L'],
                               node_size=500,
                               alpha=0.8)
        
        nx.draw_networkx_nodes(G,pos,
                               nodelist=ala_nodes,
                               node_color=node_colors['A'],
                               node_size=500,
                               alpha=0.8)
        
        nx.draw_networkx_nodes(G,pos,
                               nodelist=met_nodes,
                               node_color=node_colors['M'],
                               node_size=500,
                               alpha=0.8)

        nx.draw_networkx_nodes(G,pos,
                               nodelist=thr_nodes,
                               node_color=node_colors['T'],
                               node_size=500,
                               alpha=0.8)


        nx.draw_networkx_nodes(G,pos,
                               nodelist=lig_nodes,
                               node_color=node_colors['X'],
                               node_size=500,
                               alpha=0.8)
        
        # edges
        nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.8)        
        labels={node:node for node in G.nodes()}
        labels = {}
        nx.draw_networkx_labels(G,pos,labels,font_size=8)
        plt.savefig(figname,format='eps', dpi=300)
        #plt.show() # display       
      
    def igraph_graph(self,nodes,adjacency):
              
        if self.label_flag == "type":
            igraph_edges = [(i,nodes.index(e)) for i in range(len(nodes)) for e in adjacency[i]]
            igraph_vertices = [i for i in range(len(nodes))]
            igraph_indices = {i:nodes[i] for i in range(len(nodes))}
            igraph_labels = {i:self.igraph_color_codes[nodes[i][-1]] for i in range(len(nodes))}    #labels are given with residue type identifiers
            
            graph = ig.Graph()                 #    create igraph object
            graph.add_vertices(igraph_vertices)
            graph.add_edges(igraph_edges)
            graph.simplify(multiple = True)    #    remove double edges if there are any
            for key in igraph_labels.keys():
                graph.vs[key]["color"] = igraph_labels[key]
            
        elif self.label_flag == "assigned":
            igraph_edges = [(i,nodes.index(e)) for i in range(len(nodes)) for e in adjacency[i]]
            igraph_vertices = [i for i in range(len(nodes))]
            igraph_indices = {i:nodes[i] for i in range(len(nodes))}
            igraph_labels = {}
            igraph_color_codes_types= {}
            unlabeled = []
            for i in range(len(nodes)):
                try:
                    igraph_labels.setdefault(i,self.igraph_color_codes[nodes[i]])
                except KeyError:
                # for nodes for which labels do not exist
                    igraph_labels.setdefault(i,self.igraph_color_codes[nodes[i][-1]])   # label residue according to its type                    
            #    Get igraph object for the first graph - NOE graph
            graph = ig.Graph()                 #    create igraph object
            graph.add_vertices(igraph_vertices)
            graph.add_edges(igraph_edges)
            graph.simplify(multiple = True)    #    remove double edges if there are any
            for key in igraph_labels.keys():
                graph.vs[key]["color"] = igraph_labels[key]
        return graph,igraph_indices
       
    def community_structure(self,nodes,adjacency):
        igraph,index_igraph = self.igraph_graph(nodes,adjacency)
        community = ig.GraphBase.community_leading_eigenvector(igraph)
        community_dict = {}
        for i in range(len(community[0])):           
            community_dict.setdefault(community[0][i],[]).append(index_igraph[i])
        return community_dict    
    
    def re_order_adjecancy(self,order,old_order,old_adjacency):
        # Reorders the original adjacency list according to the new node order
        new_adjecancy = []
        for n in range(len(order)):
            adj_idx = old_order.index(order[n])
            new_adjecancy.append(old_adjacency[adj_idx])
        return new_adjecancy
               
    
    def optimal_graph_order_from_node(self,start_node,nodes,adjacency):
        # given the graph structure and the starting node (source), tries to create best order of vertices to visit
       
        edges = [(nodes[i],n) for i in range(len(nodes)) for n in adjacency[i]] # get edges from adjacency
        G = nx.Graph()
        G.add_nodes_from(nodes)
        G.add_edges_from(edges) 
        G_traverse_order = self.order_bfs_edges(nx.bfs_edges(G,start_node),start_node)
        
        while len(G_traverse_order) < len(G.nodes()):  # if not all nodes are traversed (separate subgraphs)
            leftovers = list(set(G.nodes()) - set(G_traverse_order))
            new_start_node = random.choice(leftovers)
            new_traverse_order = self.order_bfs_edges(nx.bfs_edges(G,new_start_node),new_start_node)
            G_traverse_order = G_traverse_order + new_traverse_order
        return G_traverse_order
    
    def order_bfs_edges(self,edges,start):
        order = []
        order.append(start)
        for e in edges:
            order.append(e[1])
        return order    
        
    def jaccard_sim_coefficient(self,set1,set2):
        #The Jaccard index, also known as the Jaccard similarity coefficient 
        #(originally coined coefficient by Paul Jaccard), is 
        #a statistic used for comparing the similarity and diversity of sample sets. 
        #The Jaccard coefficient measures similarity between finite sample sets, 
        #and is defined as the size of the intersection divided by the size of the union of the sample sets:
        #J (A,B) = A intersect B / A union B
        intersection, rem1, rem2 = self.extract_overlap(set1,set2)
        union = list(intersection) + list(rem1) + list(rem2)
        J = float(len(intersection))/float(len(union))
        return J
    
    def community_structure_compare(self,nodes1,adjacency1,nodes2,adjacency2):
        community_dict_1 = self.community_structure(nodes1,adjacency1)
        community_dict_2 = self.community_structure(nodes2,adjacency2)            
        #//TO DO//
        #if community graph structure will be exploited in the future
        return community_dict_1,community_dict_2
        
        
    def get_conns_dict(self,nxgraph):
        conn_dict = {}
        if len(nxgraph.edges()) > 0:
            for edge in nxgraph.edges():
                conn_dict.setdefault(edge[0],[]).append(edge[1])
                conn_dict.setdefault(edge[1],[]).append(edge[0])
            return conn_dict
        else:
            raise Exception("Contacts can not be extracted for subgraph without any edges")
         
    
    def get_patterns(self,nodes,adjacency):
        pattern_dict = {}
        for i in range(len(nodes)):
            pattern = [j[-1] for j in adjacency[i]]
            pattern = sorted(pattern)
            #print nodes[i],">>",pattern
            pattern_dict.setdefault(nodes[i],pattern)
        return pattern_dict       
    
    def compare(self,nodes1,adjacency1,nodes2,adjacency2,comparison,quantitative=False):
        #    compares neihborhood similarities for vertices of graph1 and graph2
        candidates = {}
        pattern_dict1 = self.get_patterns(nodes1,adjacency1)   
        pattern_dict2 = self.get_patterns(nodes2,adjacency2)
        if comparison == "Jaccard":
            print "Comparing vertex similarities given neighborhoods (Jaccard coefficient)... "            
            for n1 in range(len(sorted(pattern_dict1.keys()))): #,key=lambda x:int(x[:-1])))):
                for n2 in range(len(sorted(pattern_dict2.keys()))): #,key=lambda x:int(x[:-1])))):
                    # if the vertices are of the same type
                    if pattern_dict1.keys()[n1][-1] == pattern_dict2.keys()[n2][-1]:
                        sim_score = self.jaccard_sim_coefficient(pattern_dict1[pattern_dict1.keys()[n1]],pattern_dict2[pattern_dict2.keys()[n2]])
                        #print pattern_dict1.keys()[n1],">>",pattern_dict2.keys()[n2]
                        #if sim_score > 0:
                        if sim_score >= 0:
                            candidates.setdefault(pattern_dict1.keys()[n1],[]).append((pattern_dict2.keys()[n2],sim_score))   
                    else:
                        continue       
            if quantitative :
                return candidates
            else:
                candidates = self.order_candidates(candidates)
                return candidates

        #elif comparison == "Community":
            #self.community_structure_compare(nodes1,adjacency1,nodes2,adjacency2)
            #return candidates
        else:
            raise Exception("Unknown graph comparison parameter given (method compare of class Graph) ... ")
        
    def order_candidates(self,dictionary):
        #print dictionary
        for key,value in dictionary.items():
            sort_value = sorted(value,key=lambda x: x[1],reverse= True)      
            new_value = [v[0] for v in sort_value]
            dictionary[key] = new_value
        return dictionary
    
    def check_graphs_size(self,G1,G2):
        #    requirement of this algorithm is that |G1| <= |G2|; i.e. the number of vertices of G1 is smaller or equal to that of G2
        #    every node of G1 must be mapped to a node in G2
        #    "every node in G1 must be included in the correspondence"
        #    this method checks if the above requirements are fullfilled
        if len(G1.nodes())<= len(G2.nodes()):
            print 'Correct graph(s) size!'
        else:
            print sorted([n[-1] for n in G1.nodes()])
            print sorted([n[-1] for n in G2.nodes()])
            raise Exception("ERROR: |G1| must be <= |G2|!")
                            
    def check_graphs_labels(self,G1,G2):
        G1_node_labels = [node[-1] for node in G1.nodes()]
        G2_node_labels = [node[-1] for node in G2.nodes()]
        label_set1 = sorted((G1_node_labels))
        label_set2 = sorted((G2_node_labels))
        # overlap >> present in both sets
        # set1_rem >> present only in set label_set1 (G1, small graph)
        # set2_rem >> present only in set label_set2 (G2, big graph)
        overlap,set1_rem,set2_rem = self.extract_overlap(label_set1,label_set2)
        
        if len(set1_rem) == 0:
                print 'Correct labels!'
        else:
                print "Extra labels in NOE graph >> ", set1_rem
                print "Extra labels in PDB graph >> ", set2_rem
                raise Exception("ERROR: Labels mismatch!")
    
class Heuristic(Graph):
    #methods in this class provide the most optimal ordering of the graph matching priorities -> from vertices of the data graph to the vertices of the structure (pdb) graph
    def __init__(self):
        super(Graph,self).__init__()
    
    # @param
    # @retval        
    # for NOE/PDB nodes
    # given label compatibility
    # get neighboorhoods overlap
    # real value overlap assign to position in matrix
    # HungarianMaximize (matrix)   
    def hungarian_first_order(self,noe_nodes,noe_adjacency,pdb_nodes,pdb_adjacency,compatibility_dict,output_file,hops=1):
        
        fout = open(output_file,"w")
        noe_dict = self.get_conns_dict_from_adj(noe_nodes,noe_adjacency)
        pdb_dict = self.get_conns_dict_from_adj(pdb_nodes,pdb_adjacency)
        types_noe = self.sort_dict_type(noe_dict)
        types_pdb = self.sort_dict_type(pdb_dict)
        #if len(types_noe.keys())>len(types_pdb.keys()):    # crashes calculation at short distance thresholds
            #raise Exception("Difference in types of peaks and atoms >> check your input files")
        all_munkres_candidates = {}
        for type in types_noe.keys():
            candidates_dict = self.get_real_value_assign_dict(types_noe[type],types_pdb[type])
            candidate_matrix,noe_indices,pdb_indices = self.create_real_val_assign_matrix(candidates_dict)    
            Nassigned, assigned_indices = self.munkres_maximise_assignment(candidate_matrix)
            munkres_candidates = self.munkres_based_candidates(assigned_indices,candidate_matrix,noe_indices,pdb_indices)
            # write out the solution
            for key,val in munkres_candidates.items():
                fout.write("%s\t%s\n"%(key,val))
            all_munkres_candidates.update(munkres_candidates)   # add every type to the total dictionary
        fout.close()
        return all_munkres_candidates
        
    #    @params list of nodes and nested list of adjacency for 2 graph (nodes1/adjacency1 and nodes2/adjacency2)
    #    @retval    dictionary of candidate pairs where key is a node of graph1 and value is a list of tuples
    #    @retval    in values tuples:in position1 is node compatible for matching in graph2 and in position2 is a score that such matching could achieve    
    def get_real_value_assign_dict(self,conn_dict1,conn_dict2):
        # loop over all nodes in one graph        
        candidates = {}
        for peak,noeconns in conn_dict1.items():
        # for every node in a graph
            noe_neighbor_labels = [n[-1] for n in noeconns]
            for atom,pdbconns in conn_dict2.items(): # extract label compatible pdb nodes
                pdb_neighbor_labels = [p[-1] for p in pdbconns]
                common,noe_left,pdb_left = self.extract_overlap(noe_neighbor_labels,pdb_neighbor_labels)
                candidates.setdefault(peak,[]).append((atom,len(common)))
        return candidates
        
    ## Andy's heuristic method >> 10th October 2016 [original name: Craic_D]
    ## set candidate list to contain only those in dict
    # @param dic - hash table of assignment options nmr data graph vertices (keys); lists of assignment options (values)
    # @param ref - hash table of reduced assignment option for nmr data graph vertices
    # @retval - hash tabl that contains only reference assignment options for they keys in dic
    
    def update_dict(self,dic,ref):
        # IP - added the empty dictionary check
        if not bool(ref):   # in case of an emtpy reference -- return original dictionary
            return dic
                
        else:      
            for key,pars in ref.items():
                mask=np.in1d(dic[key],pars) # check if each element of reference dictionary (ref) is also element of input dictionary (dic)
                # mask is a boolean array that will select only those elements of input that are in the reference
                dic[key]=np.array(dic[key])[mask]
            return dic # 03/03/2017 Iva -- added return value
