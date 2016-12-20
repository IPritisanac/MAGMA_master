"""
MAGMA engine (c) author: Iva Pritisanac | iva.pritisanac@gmail.com;

definition and manipulation of graphs (protein structure graphs and NMR data graphs);
heuristics for ordering of the graph matching search;
munkres algorithm, variation of the Hungarian algorithm (call to the implementation given by the python package: munkres1.0.8 | https://pypi.python.org/pypi/munkres/);
algorithm for subgraph isomorphism (call to the implementation given by the high performance graph library igraph http://igraph.org/python/);
algorithm for maximal common edge subgraph (implementation is based on theoretical work from James J. McGregor;1982 | https://pdfs.semanticscholar.org/ed10/c6f788922cd5e7cb26c3d55f676979958852.pdf);
"""

import sys
import copy,collections,time
from generic import GenericMethods
import igraph as ig
import numpy as np

#import os
#import re
#from molecule import Molecule
#from myparser import Parser
#import numpy.ma as ma
#import scipy.spatial.distance as S
#import pp
#import matplotlib.pyplot as plt
#import matplotlib as mpl
#import matplotlib.cm as cm
             
class McGregor(GenericMethods):   
    # inheritance from Generic Methods
    # the class prepares the data structures for the MCES algorithm
    # the class executed MCES algorithm
    # the MCES algorithm implemented here is based on the theoretical work from McGregor J.J. 1982.

    def __init__(self,G1_node_list): 
        
        self.G1nodes = G1_node_list    
        # initialize data structures
        self.storage = [[[] for p in range(0,5)] for i in range(len(self.G1nodes)+1)]    # single storage of MEDGES/MARCS (before refinement), nodes for which matching to current node of G1 was not yet tried
        
        self.nodematch_initialpriority = [[] for node in self.G1nodes]
        self.nodematch_priority_subset = [[] for node in self.G1nodes]  # this will be dynamically updated based on the current match
        self.nodematch_all = [[] for node in self.G1nodes]
        self.allowedG2nodes = []
                
        self.currentedges = []      # edges of the currently considered G1 nodes
        self.matchededges = []      # edges of currently matched node to G1
        self.current_mapping = {}   #dynamically change --> mapping of nodes of G1 to those of G2
        
        self.edgesscore = 0             # set once MCS has been found
        self.bestedgesleft = 0          # set bestedgesleft to 0 at the beginning (once first MCS has been found, edgesleft of that MCS will become bestedgesleft)
        self.currentnode = None
        self.matchednode = None         # node matched to current node
        
        self.mcssize = 0                # size of MCS (in number of edges which are matched between the two graphs)
        self.starttime = None
        self.assignment_solutions = {} 
           
    def PrepareMcGregor(self,G1_adjacency,G2_node_list,G2_adjacency,G2_long_node_list,G2_long_adjacency,matching_priority):
        # initiate data structures
        self.G2nodes = G2_node_list       
        self.G1edges = self.EnumerateEdges(self.G1nodes,G1_adjacency)
        self.edgesleft = len(self.G1edges.keys())
        self.G2edges = self.EnumerateEdges(G2_node_list,G2_adjacency)
        self.G2edges_long = self.EnumerateEdges(G2_long_node_list,G2_long_adjacency)
        empty_medges = np.zeros(shape=(len(self.G1edges.keys()),len(self.G2edges.keys()))).astype(bool)
        empty_medges_long = np.zeros(shape=(len(self.G1edges.keys()),len(self.G2edges.keys()))).astype(bool)
        
        self.medges = self.CorrespondenceMatrix(self.G1edges,self.G2edges,empty_medges)
        self.medges_long = self.CorrespondenceMatrix(self.G1edges,self.G2edges,empty_medges_long)
        
        self.G1adjacency,self.G1indices = self.FillAdjecancyMatrix(self.G1nodes,G1_adjacency)
        self.G2adjacency,self.G2indices = self.FillAdjecancyMatrix(G2_node_list,G2_adjacency)
        self.initial_node_matchingoptions = self.SetIndexMatchPriorities(matching_priority) # initial priority dictionary can be obtained from anywhere -- this is based on the number of matching options from G2 to G1
        self.rG1indices = self.ReverseDict(self.G1indices)
        self.rG2indices = self.ReverseDict(self.G2indices)
        
        self.totG1edges = self.EdgesArrays(self.G1edges)
        self.totG2edges = self.EdgesArrays(self.G2edges)
        
    def EdgesArrays(self,edge_dict):
        edges = np.ndarray(shape=(len(edge_dict.keys()),2),dtype=int)
        for key,value in edge_dict.items():
            edges[key][0] = value[0][0]
            edges[key][1] = value[0][1]
        return edges            
        
    def ReverseDict(self,init_dict):
        # given dictionary with single key and value, returns reversed dictionary (single value: key)
        return {value:key for key,value in init_dict.items()}
                               
    def FillAdjecancyMatrix(self,nodes,adjacency):
        node_indexing = {}        
        adjacency_matrix = np.zeros(shape=(len(nodes),len(nodes))) 
        for node in range(len(nodes)):            # gets to node index
            node_indexing.setdefault(nodes[node],node)
            for neighbour in adjacency[node]:     # gets to node neighbours
                neighbour_index = nodes.index(neighbour)
                adjacency_matrix[node,neighbour_index] = 1
        return adjacency_matrix,node_indexing
    
    def SetIndexMatchPriorities(self,match_prior):
        node_match_options = {}        
        for G1,G2s in match_prior.items():
            G2indices = [self.G2indices[g] for g in G2s]
            node_match_options.setdefault(self.G1indices[G1],G2indices)
        return node_match_options

    def EnumerateNodes(self,nodes):
        enumerated_nodes = {}
        for n in range(len(nodes)):    #n starts from 0
            enumerated_nodes.setdefault(n,nodes[n])
        return enumerated_nodes
    
    def EnumerateEdges(self,nodes,adjacency):
        edges={}
        enumerated_edges = {}
        for i in range(len(adjacency)):
            node=nodes[i]    #node currently considered
            typei=node[-1:]  #type of node i
            for j in adjacency[i]:
                indexj=nodes.index(j)
                key=(i,indexj)      #indexes of each of the nodes in nodes list
                typej=j[-1:]        
                tup=(typei,typej)   #types of each of the nodes
                listtup=sorted(tup) #sort types tuple
                edges.setdefault(key,tuple(listtup))   #edges type tuples are sorted -- easier comparison to edges of another graph (no double counting will be needed)

        cnte=-1    #counter of edges from adjacency list, ensure that counting starts from 0
        unique_edges={}
        exclude=[]    #list of reverse keys in edges (the same edge) to be excluded if exists elsewhere in the dictionary
        for key,value in edges.items():
            exclude.append((key[1],key[0]))    #append inverse of the key, so that same edge is not double counted
            if key in exclude:
                continue
            else:
                cnte+=1
                unique_edges.setdefault(key,value)
                enumerated_edges.setdefault(cnte,(key,value))
        return enumerated_edges
    
    def Neighbourhood(self,nodes,adjacency):
        neighbourhood = {}
        for i in range(len(adjacency)):
            pattern=[]
            if len(adjacency[i])>0:
                for j in adjacency[i]:
                    pattern.append(j[-1:])
                    sort=sorted(pattern)   #sorting allows for pattern comparison between graphs
                neighbourhood.setdefault(nodes[i],tuple(sort))
            else:
                print 'No pattern for ', i
                neighbourhood.setdefault(nodes[i],tuple())
        return neighbourhood
    
    def CorrespondenceMatrix(self,edges1,edges2,empty_matrix):
        pairs=[]
        for key1,value1 in edges1.items():
            for key2,value2 in edges2.items():    #values are sorted (i.e. edge (I,V)==(V,I) would be (I,V))
                if value1[1]==value2[1]:
                    pairs.append((key1,key2))
        for pair in pairs:
            empty_matrix[pair[0],pair[1]]=1             
        return empty_matrix

    def GetEdges(self):
        #    information from initial edges-matrix could be used, without the need to check types!
        #    emptied every time the function is called
        #    get all edges of G1 node
        self.currentedges = np.where(self.totG1edges == self.G1node)[0]       
        self.matchededges = np.where(self.totG2edges == self.matchednode)[0]        
        
    def UpdateMEDGES(self):
        
        mask1 = np.zeros_like(self.medges).astype(bool) 
        mask2 = np.ones_like(self.medges).astype(bool)
        mask1[self.currentedges]=True
        mask2[:,self.matchededges]=False
        self.medges[np.logical_and(mask1,mask2)]=False

    def SetPriorities(self):        
        # based on the tentative mapping (self.G1 node >> self.matchednode), set priority lists for other nodes
        #for each node k of G, such that k > i and k is adjacent to node i do
        #priority subset [K] := priority subset [l] n {l :l is adjacent to j in G,).

        adjacent_i = np.nonzero(self.G1adjacency[self.G1node])[0]    # indices of i adjacent nodes
        adjacent_j = np.nonzero(self.G2adjacency[self.matchednode])[0]    # indices of j adjacent nodes
        
        for k in adjacent_i:
            if k > self.G1node:
                #for k in range(len(self.G1nodes)):  
                #    if k > self.G1node and k in adjacent_i:   # if k is subsequent to current node in G1 in order in which nodes are matched
                                
                priority_subset = [l for l in adjacent_j if l not in self.current_mapping.values() and l in self.initial_node_matchingoptions[k]]     # and not yet assigned
                                
                a_multiset = collections.Counter(self.nodematch_priority_subset[k])
                b_multiset = collections.Counter(priority_subset)
                overlap = list((a_multiset & b_multiset).elements())       
                self.nodematch_priority_subset[k] = self.ScoreAdjust(k,overlap)
                                    
                #self.nodematch_priority_subset[k] = bestoption                
                
                all_others = [node for node in self.nodematch_initialpriority[k] if node not in overlap]                    
                self.nodematch_all[k] = self.ScoreAdjust(k,all_others) 
        
    ## This function stores priority matching
    # @retval        
    def StorePriorities(self):      
        storeindx = self.G1node + 1 
        self.storage[storeindx][0] = copy.deepcopy(self.nodematch_priority_subset)
        self.storage[storeindx][1] = copy.deepcopy(self.nodematch_all)
        
    def StoreMedgesEdgesleft(self):
        # when last node is reached forward storing cannot be done
        storeindx = self.G1node + 1
        self.storage[storeindx][3] = copy.deepcopy(self.medges)  # store initial medges/edgesleft to workspace connected with node 0
        self.storage[storeindx][4] = copy.deepcopy(self.edgesleft)  

    def StoreInitial(self):
        # important for when backtrack to starting node happens --> values before any assignment to node 0
        self.storage[0][0] = copy.deepcopy(self.nodematch_priority_subset)
        self.storage[0][1] = copy.deepcopy(self.nodematch_all)
        self.storage[0][3] = copy.deepcopy(self.medges) # store initial medges/edgesleft to workspace connected with node 0
        self.storage[0][4] = copy.deepcopy(self.edgesleft)

    def StoreTried(self):
        self.storage[self.G1node][2].append(copy.deepcopy(self.matchednode))   # mark node xi as tried for node i
              
    def InitialG2PriorityMatchList(self):     
        #    set initial node order (order in which nodes of G1 will be picked)
        #    initial options for matching of nodes in G1 --> critical is that indexing is order same way as the treesearch (treesearch == G1 nodes)
        for n1,priority_list in self.initial_node_matchingoptions.items():
            self.nodematch_priority_subset[n1] = priority_list
            self.nodematch_initialpriority[n1] = priority_list
            self.nodematch_all[n1] = priority_list
                                        
    def MatchOptions(self):  
        #    dynamically update possible matches considering current mapping and previously tried ones
        #    check the storage of currently tried match options ( in pseudocode (if there are any UNTRIED options to which G1 may correspond to))       
        #print self.G1node
        #print len(self.nodematch_priority_subset)
        priorityG2options = self.nodematch_priority_subset[self.G1node]   # all priority G2 options
        otherG2options = self.nodematch_all[self.G1node]    # leftover G1 options
        triedG2 = self.storage[self.G1node][2]   # all currently untried G2 options    
        
        allowedG2nodes = [node for node in priorityG2options if node not in self.current_mapping.values() and node not in triedG2]    # G2 nodes to which i of G1 may correspond to
        bestoptions = self.ScoreAdjust(self.G1node,allowedG2nodes)
        self.allowedG2nodes = bestoptions
        
        if len(allowedG2nodes) == 0:
            allowedG2nodes = [node for node in otherG2options if node not in self.current_mapping.values() and node not in triedG2]    # G2 nodes to which i of G1 may correspond to       
            bestoptions = self.ScoreAdjust(self.G1node,allowedG2nodes)
            self.allowedG2nodes = bestoptions
                     
    def ScoreAdjust(self,G1node,options):
        sorted_options = []
        for G2node in options:
            indexG2 = self.nodematch_initialpriority[G1node].index(G2node)
            sorted_options.append((G2node,indexG2))
        sorted_options = sorted(sorted_options,key=lambda x: x[1])
        sorted_options = [sorted_options[i][0] for i in range(len(sorted_options))]
        # print >>sorted_options
        return sorted_options
    
    def RePrioritize(self,assigned,init_matchingoptions):

        new_priorities = {}
        for peak,atoms in init_matchingoptions.items():
            old_prior = [a for a in atoms if a not in assigned[peak]]
            new_prior = assigned[peak]+old_prior
            new_priorities.setdefault(peak,new_prior)
        print "Reordered match order >> ",new_priorities
        return new_priorities                   
      
    def CollectBestAssignments(self,current_solution):
        self.assignment_solutions.setdefault(self.edgesleft,{}) # sort assignments according to their score
        for key,value in current_solution.items():
            self.assignment_solutions[self.edgesleft].setdefault(self.rG1indices[key],[])
            if self.rG2indices[value] not in self.assignment_solutions[self.edgesleft][self.rG1indices[key]]: # if this is a new assignment
                self.assignment_solutions[self.edgesleft][self.rG1indices[key]].append(self.rG2indices[value])
    
    def GetFinalAssignments(self):
        final_assignments = {}
        for score in sorted(self.assignment_solutions.keys(),reverse=True):
            for peak,atoms in self.assignment_solutions[score].items():
                final_assignments.setdefault(peak,[])
                for atom in atoms:
                    if atom not in final_assignments[peak]:
                        final_assignments[peak].append(atom)
        return final_assignments   
      
    def McGregorLoop(self,outfile,runall=True,n_mces=None,time_check=False,maximum_time=None,thresh=True):
        self.InitialG2PriorityMatchList()   # set initial priority matching list              
        self.G1node = 0  #    index of start node in G1.nodes list
        self.StoreInitial()
        #store initial priorities, MEDGES, EDGESLEFT (no refinements made) -- stored in the workspace of node 0 of G1 graph  
        self.currentnode = self.G1nodes[self.G1node]     # start at the beginning of the list
        self.starttime = time.time()
        self.output = open(outfile,"w") # open the output file where resulting assignments will be stored
        iter_cnt = 0    # count how many exact MCES iterations performed
        mces_cnt = 0    # count how many MCES have been found        
        while True:             
            self.MatchOptions()     # sets the lists self.untriedG2nodes and self.allowedG2nodes for possible matches to current node 
            iter_cnt+=1	#	count an iteration of MCES algorithm       
            #########
            # Test terminating condition for distance thresholding
            # if maximum number of iterations or maximum time is reached >> exit
            # max time is set to 24 hours
            #if thresh and int(time.time()-self.starttime) > 86400:
            #if thresh and iter_cnt > len(self.G1nodes)*1000 or int(time.time()-self.starttime) > 86400:
            #print "N of iterations or time exceeded maximum ... exiting"
            #sys.exit(100)
            #########
            if  len(self.allowedG2nodes)!=0:
                self.matchednode = self.allowedG2nodes[0]   # take the first (best) option,from allowed nodes, pick the node from G2, which would be the best option for current G1 node                
                self.current_mapping[self.G1node] = self.matchednode   # keeps track of current mappings                
                self.StoreTried()
                self.GetEdges()   #    get to edges of current node in G1 and its matched node in G2
                # before refinement of MEDGES on the basis of this correspondence, create temporary copy of MEDGES which will either be accepted or rejected
                medges_tmp = copy.deepcopy(self.medges)
                # REFINE MEDGES on the basis of this tentative correspondence
                self.UpdateMEDGES()
                self.edgesleft = np.sum(np.any(self.medges,axis=1))               
                self.SetPriorities()
                self.StorePriorities()  # store these priorities in the workspace of i + 1 node of G1                                
                # in pseudocode -- > if there are untried! nodes in G2 to which node of G1 may(not matched to others already) correspond to, then Xi := one of these nodes, mark G2 as tried for i
                if self.edgesleft > self.bestedgesleft or (runall and self.edgesleft >= self.bestedgesleft):
                    if self.G1node == self.G1nodes.index(self.G1nodes[-1]): # if we came to the end of the node list 
                        print "found MCES of size >> ",self.edgesleft,"in iter >> ",iter_cnt
                        mces_cnt+=1  
                        current = copy.deepcopy(self.current_mapping)
                        self.CollectBestAssignments(current)
                        self.bestedgesleft = copy.deepcopy(self.edgesleft)    # dynamically assign bestedges score to edgesleft score every time MCS is found 

                        # look up leftover assignment options 
                        # commented out on 26/10/16 - Iva
                        """
                        nodeleft=np.setdiff1d(self.initial_node_matchingoptions[self.G1node],self.current_mapping.values())
                        if len(nodeleft) > 0: # all of its assignment options
                            #entering here means that there are "STILL SOME OPTIONS!",nodeleft
                            currmap = copy.deepcopy(self.current_mapping)
                            for n in nodeleft:
                                currmap[self.G1node] = n
                                endmedges = copy.deepcopy(self.storage[self.G1node-1][3])
                                endedgesleft = copy.deepcopy(self.storage[self.G1node-1][4])
                                mask1 = np.zeros_like(endmedges).astype(bool) 
                                mask2 = np.ones_like(endmedges).astype(bool)
                                edges1 = np.where(self.totG1edges == self.G1node)[0]    
                                edges2 = np.where(self.totG2edges == n)[0]
                                mask1[edges1]=True
                                mask2[:,edges2]=False
                                endmedges[np.logical_and(mask1,mask2)]=False
                                endedgesleft = np.sum(np.any(endmedges,axis=1))
                                if endedgesleft > self.bestedgesleft or (runall and endedgesleft >= self.bestedgesleft):
                                    #print "found MCES of size >> ",endedgesleft,"in iter >> ",iter_cnt
                                    mces_cnt+=1
                                    self.bestedgesleft = endedgesleft
                                    self.edgesleft = endedgesleft
                                    #print endedgesleft," Found another MCES!"
                                    self.CollectBestAssignments(currmap)
                                    self.output.write('%s\n'%(endedgesleft))
                                    for key,value in currmap.items():
                                        self.output.write('%s\t%s\n'%(self.rG1indices[key],self.rG2indices[value]))
                                    self.output.write('\n')
                                    self.output.flush()
                                    elapsed=time.time()-self.starttime                                    
                                    if (n_mces and mces_cnt >= n_mces) or (time_check and elapsed>maximum_time):
                                        self.output.close()
                                        final_assign_solutions = self.GetFinalAssignments()
                                        return self.bestedgesleft,final_assign_solutions
                                        """
                        # indicates that MCS has been found --> this will influence the treeorder appending
                        #self.bestedgesleft = copy.deepcopy(self.edgesleft)    # dynamically assign bestedges score to edgesleft score every time MCS is found 
                        timeittook = time.time() - self.starttime
                        self.output.write('%s\n'%(timeittook))
                        self.output.write('%s\n'%(self.edgesleft))
                                               
                        for key,value in self.current_mapping.items():                          
                            self.output.write('%s\t%s\n'%(self.rG1indices[key],self.rG2indices[value]))
                        self.output.write('\n')
                        self.output.flush()
                        self.StoreMedgesEdgesleft()
                        elapsed=time.time()-self.starttime                        
                        if (n_mces and mces_cnt >= n_mces) or (time_check and elapsed>maximum_time):
                            self.output.close()
                            final_assign_solutions = self.GetFinalAssignments()
                            return self.bestedgesleft,final_assign_solutions
                    else:  
                        self.StoreMedgesEdgesleft()
                        self.G1node= self.G1node + 1
                        self.storage[self.G1node][2]=[] # mark all nodes as untried
                else:
                    # bring back unchanges medges for next iteration / G2 option checking
                    self.medges = medges_tmp
            else:       
                # if all matching options for current node have been tried
                # move a node up in a tree search       
                #print "backtracking", self.G1node        
                if self.G1node in self.current_mapping.keys():
                    self.current_mapping.pop(self.G1node) 
                self.G1node = self.G1node-1              
                self.medges = copy.deepcopy(self.storage[self.G1node][3])
                self.edgesleft = copy.deepcopy(self.storage[self.G1node][4])
                self.nodematch_all = copy.deepcopy(self.storage[self.G1node][1])
                self.nodematch_priority_subset = copy.deepcopy(self.storage[self.G1node][0])
                #print 'Edgesleft ', self.edgesleft
            # Iva - 20/12/2016 -- added the condition self.allowedG2nodes==0
            if (self.G1node==-1) and (self.allowedG2nodes==0):
                break
        self.output.close()
        final_assign_solutions = self.GetFinalAssignments()
        return self.bestedgesleft,final_assign_solutions                
