"""
__authors__: Iva Pritisanac and Andy Baldwin

implementation of McGregor algorithm for maximal common edge subgraph in

        > in c (class McGregor; AJB)
        > in python (class MCES_PY; IP)

implementation based on theoretical work from James J. McGregor;1982 | 
https://pdfs.semanticscholar.org/ed10/c6f788922cd5e7cb26c3d55f676979958852.pdf);
"""

import copy,collections,time
from generic import GenericMethods
import numpy as np
import os,sys

class McGregor(GenericMethods):   # AJB Oct 2016 version of McGregor class for C version of the algorithm on the basis of python class MCES_PY

    def __init__(self,G1_node_list): 
        
        self.G1nodes = G1_node_list    
        
    def prepare_mcgregor(self,G1_adjacency,G2_node_list,G2_adjacency,G2_long_node_list,G2_long_adjacency,matching_priority,outdir):
        
        self.outdir=outdir
        # IP 03/03/2017 not necessary -- checked for by Magma class
        #if(os.path.exists(self.outdir)==0):
            #os.system('mkdir '+self.outdir)
            
        # IP 03/03/2017 -- adjusted to make system independent
        if not os.path.exists(self.outdir+os.sep+"core"): # if the output results directory does not exist on this path
            os.makedirs(self.outdir+os.sep+"core")   # create the directory
        
        outy=open(self.outdir+os.sep+"core"+os.sep+"mcesCore.init",'w') #self.outdir+'/core/mcesCore.init','w')
        outy.write('G1_nodes\n')
        for i in range(len(self.G1nodes)):
            outy.write('%s :\t' % (self.G1nodes[i]))
            for j in range(len(G1_adjacency[i])):
                outy.write('%s\t' % (G1_adjacency[i][j]))
            outy.write('\n')
        outy.write('G2_nodes\n')
        for i in range(len(G2_node_list)):
            outy.write('%s :\t' % (G2_node_list[i]))
            for j in range(len(G2_adjacency[i])):
                outy.write('%s\t' % (G2_adjacency[i][j]))
            outy.write('\n')
        outy.write('MatchPriorities\n')
        for i in range(len(self.G1nodes)):
            outy.write('%s :\t' % (self.G1nodes[i]))
            for j in range(len(matching_priority[self.G1nodes[i]])):
                outy.write('%s\t' % (matching_priority[self.G1nodes[i]][j]))
            outy.write('\n')
        outy.close()

        #for item,values in matching_priority.items():
        #print item,':',values

    ## create a new vertex matching priority list for every vertex in the dictionary on the basis of assigned
    # @param assigned - updated hash table where every vertex (key) has with it associated list of assignment options (values)
    # @param init_matchingoptions - original hash table where every vertex (key) has with it associated list of assignment options (values)
    # @retval new_priorities - an updated hash table of assignment options; or in cases of errors in updating, the original init_matchingoptions 
    def mc_reprioritize(self,assigned,init_matchingoptions):
        
        if not bool(assigned):  # if the entering assigned dictionary is empty
            print "An empty assignment dictionary entered in mc_reprioritize assigned; returning the original dictionary of assignment options"
            return init_matchingoptions # return the original candidate matching options
        elif len(assigned.keys())!=len(init_matchingoptions.keys()): # if the entering dictionary does not have an assignment option for all vertices
            print "The assignment dictionary entering in mc_reprioritize does not match the original; returning the original dictionary of assignment options"
            return init_matchingoptions
        else:
            try:
                new_priorities = {}
                for peak,atoms in init_matchingoptions.items():
                    old_prior = [a for a in atoms if a not in assigned[peak]]
                    new_prior = assigned[peak]+old_prior    # set the beginning of the list to the entering assigned values, append to the end the original values
                    new_priorities.setdefault(peak,new_prior)
                    #print "Reordered match order >> " #,new_priorities
                #for key,value in new_priorities.items():
                    #print key,':',value
                return new_priorities
            except ValueError:  # in case keys of the two dictionaries don't match
                print "mc_reprioritize caught ValueError, returning the original dictionary of assignment options"
                return init_matchingoptions
        
        
    def mcgregor_loop(self,outfile,runall=True,n_mces=None,time_check=False,maximum_time=None,thresh=True,parflg=0): # AJB Oct 2016 
            
        try: #this is the number of mces-es to acquire. 0 means get all.
            n_mcesSet=int(n_mces)
        except:
            n_mcesSet=0

        if(runall==True): #save all above a threshold size or just the first?
            runSet=1
        else:
            runSet=0
            
        if(outfile==False):
            runOut=str(0)
        else:
            runOut=outfile
            
        try:
            maxtimeSet=int(maximum_time)
        except:
            maxtimeSet=0

        # IP 04/03/2017
        # assumes that the input text file is located inside dir MAGMA/ 
        # assumes that bin is located at path/to/MAGMA/bin        
        #cpath = ".."+os.sep+".."+os.sep+"src"+os.sep+"mcesCore"
        cpath = os.path.abspath("bin"+os.sep+"mcesCore") # get the absolute path to the bin directory and the core of the C code for the mces algorithm
        if(parflg==0):  
            runline = cpath+" "+self.outdir+os.sep+"core"+os.sep+"mcesCore.init "+runOut+" "+self.outdir+os.sep+"core"+os.sep+"final.out "+str(n_mcesSet)+" "+str(runSet)+" "+str(maxtimeSet)+" "+str(parflg) 
            #runline='../../src/mcesCore '+self.outdir+'/core/mcesCore.init '+runOut+' '+self.outdir+'/core/final.out '+str(n_mcesSet)+' '+str(runSet)+' '+str(maxtimeSet)+' '+str(parflg)
        else:
            runline = cpath+" "+self.outdir+os.sep+"core"+os.sep+"mcesCore.init "+runOut+" "+self.outdir+os.sep+"core"+os.sep+"final.out "+str(n_mcesSet)+" "+str(runSet)+" "+str(maxtimeSet)+" "+str(0)
            #runline='mpiexec -np 2 ../../src/mcesCore '+self.outdir+'/core/mcesCore.init '+runOut+' '+self.outdir+'/core/final.out '+str(n_mcesSet)+' '+str(runSet)+' '+str(maxtimeSet)+' '+str(parflg)
            #runline='../../src/mcesCore '+self.outdir+'/core/mcesCore.init '+runOut+' '+self.outdir+'/core/final.out '+str(n_mcesSet)+' '+str(runSet)+' '+str(maxtimeSet)+' '+str(0)

        sys.stdout.flush()

        os.system(runline) # does not work on windows
        
        if(runSet==1):
            print runline
            
        #read in output file
        if(os.path.exists(self.outdir+os.sep+"core"+os.sep+"final.out")==1): #'/core/final.out')==1):
            inny=open(self.outdir+os.sep+"core"+os.sep+"final.out")
            final_assignments={}
            for line in inny.readlines():
                test=line.split()
                if(len(test)==1):
                    edgesleft=int(test[0])
                else:
                    key=test[0].split(':')[0]
                    ass=[]
                    for j in range(len(test)-1):
                        ass.append(test[j+1])
                    final_assignments[key]=ass
            os.system('rm -rf '+self.outdir+os.sep+'core')
            return edgesleft,final_assignments
        else:
            print 'No output'
            return 0,{}
        
        #if(edgesleft==0):
        #    print 'Problem: no score. Aborting'
        #    sys.exit(1)             

# inheritance from Generic Methods
# the class prepares the data structures for the MCES algorithm
# the class executed MCES algorithm
# the MCES algorithm implemented here is based on the theoretical work from McGregor J.J. 1982.
class MCES_PY(GenericMethods):   # IP python version
    
    def __init__(self,g1_node_list): 
        
        self.g1_nodes = g1_node_list    
        # initialize data structures
        self.storage = [[[] for p in range(0,5)] for i in range(len(self.g1_nodes)+1)]    # single storage of MEDGES/MARCS (before refinement), nodes for which matching to current node of G1 was not yet tried
        
        self.nodematch_initialpriority = [[] for node in self.g1_nodes]
        self.nodematch_priority_subset = [[] for node in self.g1_nodes]  # this will be dynamically updated based on the current match
        self.nodematch_all = [[] for node in self.g1_nodes]
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
           
    def prepare_mces_py(self,G1_adjacency,G2_node_list,G2_adjacency,G2_long_node_list,G2_long_adjacency,matching_priority):
        # initiate data structures
        self.G2nodes = G2_node_list       
        self.G1edges = self.enumerate_edges(self.g1_nodes,G1_adjacency)
        self.edgesleft = len(self.G1edges.keys())
        self.G2edges = self.enumerate_edges(G2_node_list,G2_adjacency)
        self.G2edges_long = self.enumerate_edges(G2_long_node_list,G2_long_adjacency)
        empty_medges = np.zeros(shape=(len(self.G1edges.keys()),len(self.G2edges.keys()))).astype(bool)
        empty_medges_long = np.zeros(shape=(len(self.G1edges.keys()),len(self.G2edges.keys()))).astype(bool)
        
        self.medges = self.correspondence_matrix(self.G1edges,self.G2edges,empty_medges)
        self.medges_long = self.correspondence_matrix(self.G1edges,self.G2edges,empty_medges_long)
        
        self.G1adjacency,self.G1indices = self.fill_adjacency_matrix(self.g1_nodes,G1_adjacency)
        self.G2adjacency,self.G2indices = self.fill_adjacency_matrix(G2_node_list,G2_adjacency)
        self.initial_node_matchingoptions = self.set_index_match_priorities(matching_priority) # initial priority dictionary can be obtained from anywhere -- this is based on the number of matching options from G2 to G1
        self.rG1indices = self.reverse_dict(self.G1indices)
        self.rG2indices = self.reverse_dict(self.G2indices)
        
        self.totG1edges = self.edges_array(self.G1edges)
        self.totG2edges = self.edges_array(self.G2edges)
        
    def edges_array(self,edge_dict):
        edges = np.ndarray(shape=(len(edge_dict.keys()),2),dtype=int)
        for key,value in edge_dict.items():
            edges[key][0] = value[0][0]
            edges[key][1] = value[0][1]
        return edges            
        
    def reverse_dict(self,init_dict):
        # given dictionary with single key and value, returns reversed dictionary (single value: key)
        return {value:key for key,value in init_dict.items()}
                               
    def fill_adjacency_matrix(self,nodes,adjacency):
        node_indexing = {}        
        adjacency_matrix = np.zeros(shape=(len(nodes),len(nodes))) 
        for node in range(len(nodes)):            # gets to node index
            node_indexing.setdefault(nodes[node],node)
            for neighbour in adjacency[node]:     # gets to node neighbours
                neighbour_index = nodes.index(neighbour)
                adjacency_matrix[node,neighbour_index] = 1
        return adjacency_matrix,node_indexing
    
    def set_index_match_priorities(self,match_prior):
        node_match_options = {}        
        for G1,G2s in match_prior.items():
            G2indices = [self.G2indices[g] for g in G2s]
            node_match_options.setdefault(self.G1indices[G1],G2indices)
        return node_match_options

    def enumerate_nodes(self,nodes):
        enumerated_nodes = {}
        for n in range(len(nodes)):    #n starts from 0
            enumerated_nodes.setdefault(n,nodes[n])
        return enumerated_nodes
    
    def enumerate_edges(self,nodes,adjacency):
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
    
    """
    //not used//
    """
    def define_neighbourhood(self,nodes,adjacency):
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
    
    def correspondence_matrix(self,edges1,edges2,empty_matrix):
        pairs=[]
        for key1,value1 in edges1.items():
            for key2,value2 in edges2.items():    #values are sorted (i.e. edge (I,V)==(V,I) would be (I,V))
                if value1[1]==value2[1]:
                    pairs.append((key1,key2))
        for pair in pairs:
            empty_matrix[pair[0],pair[1]]=1             
        return empty_matrix

    def get_edges(self):
        #    information from initial edges-matrix could be used, without the need to check types!
        #    emptied every time the function is called
        #    get all edges of G1 node
        self.currentedges = np.where(self.totG1edges == self.G1node)[0]       
        self.matchededges = np.where(self.totG2edges == self.matchednode)[0]        
        
    def update_medges(self):
        
        mask1 = np.zeros_like(self.medges).astype(bool) 
        mask2 = np.ones_like(self.medges).astype(bool)
        mask1[self.currentedges]=True
        mask2[:,self.matchededges]=False
        self.medges[np.logical_and(mask1,mask2)]=False

    def set_priorities(self):        
        # based on the tentative mapping (self.G1 node >> self.matchednode), set priority lists for other nodes
        #for each node k of G, such that k > i and k is adjacent to node i do
        #priority subset [K] := priority subset [l] n {l :l is adjacent to j in G,).

        adjacent_i = np.nonzero(self.G1adjacency[self.G1node])[0]    # indices of i adjacent nodes
        adjacent_j = np.nonzero(self.G2adjacency[self.matchednode])[0]    # indices of j adjacent nodes
        
        for k in adjacent_i:
            if k > self.G1node:
                #for k in range(len(self.g1_nodes)):  
                #    if k > self.G1node and k in adjacent_i:   # if k is subsequent to current node in G1 in order in which nodes are matched
                                
                priority_subset = [l for l in adjacent_j if l not in self.current_mapping.values() and l in self.initial_node_matchingoptions[k]]     # and not yet assigned
                                
                a_multiset = collections.Counter(self.nodematch_priority_subset[k])
                b_multiset = collections.Counter(priority_subset)
                overlap = list((a_multiset & b_multiset).elements())       
                self.nodematch_priority_subset[k] = self.score_adjust(k,overlap)
                                    
                #self.nodematch_priority_subset[k] = bestoption                
                
                all_others = [node for node in self.nodematch_initialpriority[k] if node not in overlap]                    
                self.nodematch_all[k] = self.score_adjust(k,all_others) 
        
    ## This function stores priority matching
    # @retval        
    def store_priorities(self):      
        storeindx = self.G1node + 1 
        self.storage[storeindx][0] = copy.deepcopy(self.nodematch_priority_subset)
        self.storage[storeindx][1] = copy.deepcopy(self.nodematch_all)
        
    def store_medges_edgesleft(self):
        # when last node is reached forward storing cannot be done
        storeindx = self.G1node + 1
        self.storage[storeindx][3] = copy.deepcopy(self.medges)  # store initial medges/edgesleft to workspace connected with node 0
        self.storage[storeindx][4] = copy.deepcopy(self.edgesleft)  

    def store_initial(self):
        # important for when backtrack to starting node happens --> values before any assignment to node 0
        self.storage[0][0] = copy.deepcopy(self.nodematch_priority_subset)
        self.storage[0][1] = copy.deepcopy(self.nodematch_all)
        self.storage[0][3] = copy.deepcopy(self.medges) # store initial medges/edgesleft to workspace connected with node 0
        self.storage[0][4] = copy.deepcopy(self.edgesleft)

    def store_tried(self):
        self.storage[self.G1node][2].append(copy.deepcopy(self.matchednode))   # mark node xi as tried for node i
              
    def initial_g2_priority_match_list(self):     
        #    set initial node order (order in which nodes of G1 will be picked)
        #    initial options for matching of nodes in G1 --> critical is that indexing is order same way as the treesearch (treesearch == G1 nodes)
        for n1,priority_list in self.initial_node_matchingoptions.items():
            self.nodematch_priority_subset[n1] = priority_list
            self.nodematch_initialpriority[n1] = priority_list
            self.nodematch_all[n1] = priority_list
                                        
    def match_options(self):  
        #    dynamically update possible matches considering current mapping and previously tried ones
        #    check the storage of currently tried match options ( in pseudocode (if there are any UNTRIED options to which G1 may correspond to))       
        #print self.G1node
        #print len(self.nodematch_priority_subset)
        priorityG2options = self.nodematch_priority_subset[self.G1node]   # all priority G2 options
        otherG2options = self.nodematch_all[self.G1node]    # leftover G1 options
        triedG2 = self.storage[self.G1node][2]   # all currently untried G2 options    
        
        allowedG2nodes = [node for node in priorityG2options if node not in self.current_mapping.values() and node not in triedG2]    # G2 nodes to which i of G1 may correspond to
        bestoptions = self.score_adjust(self.G1node,allowedG2nodes)
        self.allowedG2nodes = bestoptions
        
        if len(allowedG2nodes) == 0:
            allowedG2nodes = [node for node in otherG2options if node not in self.current_mapping.values() and node not in triedG2]    # G2 nodes to which i of G1 may correspond to       
            bestoptions = self.score_adjust(self.G1node,allowedG2nodes)
            self.allowedG2nodes = bestoptions
                     
    def score_adjust(self,G1node,options):
        sorted_options = []
        for G2node in options:
            indexG2 = self.nodematch_initialpriority[G1node].index(G2node)
            sorted_options.append((G2node,indexG2))
        sorted_options = sorted(sorted_options,key=lambda x: x[1])
        sorted_options = [sorted_options[i][0] for i in range(len(sorted_options))]
        # print >>sorted_options
        return sorted_options
    
    """
    This method bugs occassionally -- check
    """
    def re_prioritize(self,assigned,init_matchingoptions):
        
        if not bool(assigned):  # if the entering assigned dictionary is empty
            print "An empty assignment dictionary entered in re_prioritize assigned; returning the original dictionary of assignment options"
            return init_matchingoptions # return the original candidate matching options

        elif len(assigned.keys())!=len(init_matchingoptions.keys()): # if the entering dictionary does not have an assignment option for all vertices
            print "The assignment dictionary entering in re_prioritize does not match the original; returning the original dictionary of assignment options"
            return init_matchingoptions
        else:
            try:
                new_priorities = {}
                for peak,atoms in init_matchingoptions.items():
                    old_prior = [a for a in atoms if a not in assigned[peak]]
                    new_prior = assigned[peak]+old_prior    # set the beginning of the list to the entering assigned values, append to the end the original values
                    new_priorities.setdefault(peak,new_prior)
                    #print "Reordered match order >> " #,new_priorities
                #for key,value in new_priorities.items():
                    #print key,':',value
                return new_priorities
            except ValueError:  # in case keys of the two dictionaries don't match
                print "re_prioritize caught ValueError, returning the original dictionary of assignment options"
                return init_matchingoptions                 
      
    def collect_best_assignments(self,current_solution):
        self.assignment_solutions.setdefault(self.edgesleft,{}) # sort assignments according to their score
        for key,value in current_solution.items():
            self.assignment_solutions[self.edgesleft].setdefault(self.rG1indices[key],[])
            if self.rG2indices[value] not in self.assignment_solutions[self.edgesleft][self.rG1indices[key]]: # if this is a new assignment
                self.assignment_solutions[self.edgesleft][self.rG1indices[key]].append(self.rG2indices[value])
    
    def get_final_assignments(self):
        final_assignments = {}
        for score in sorted(self.assignment_solutions.keys(),reverse=True):
            for peak,atoms in self.assignment_solutions[score].items():
                final_assignments.setdefault(peak,[])
                for atom in atoms:
                    if atom not in final_assignments[peak]:
                        final_assignments[peak].append(atom)
        return final_assignments   
      
    def mcgregor_main(self,outfile,runall=True,n_mces=None,time_check=False,maximum_time=None,thresh=True):
        
        self.initial_g2_priority_match_list()   # set initial priority matching list              
        self.G1node = 0  #    index of start node in G1.nodes list
        self.store_initial()
        #store initial priorities, MEDGES, EDGESLEFT (no refinements made) -- stored in the workspace of node 0 of G1 graph  
        self.currentnode = self.g1_nodes[self.G1node]     # start at the beginning of the list
        self.starttime = time.time()
        self.output = open(outfile,"w") # open the output file where resulting assignments will be stored
        iter_cnt = 0    # count how many exact MCES iterations performed
        mces_cnt = 0    # count how many MCES have been found        
        while True:             
            self.match_options()     # sets the lists self.untriedG2nodes and self.allowedG2nodes for possible matches to current node 
            iter_cnt+=1	#	count an iteration of MCES algorithm
            #########
            if  len(self.allowedG2nodes)!=0:
                self.matchednode = self.allowedG2nodes[0]   # take the first (best) option,from allowed nodes, pick the node from G2, which would be the best option for current G1 node                
                self.current_mapping[self.G1node] = self.matchednode   # keeps track of current mappings                
                self.store_tried()
                self.get_edges()   #    get to edges of current node in G1 and its matched node in G2
                # before refinement of MEDGES on the basis of this correspondence, create temporary copy of MEDGES which will either be accepted or rejected
                medges_tmp = copy.deepcopy(self.medges)
                # REFINE MEDGES on the basis of this tentative correspondence
                self.update_medges()
                self.edgesleft = np.sum(np.any(self.medges,axis=1))               
                self.set_priorities()
                self.store_priorities()  # store these priorities in the workspace of i + 1 node of G1                                
                ## if there are untried! nodes in G2 to which node of G1 may(not matched to others already) correspond to, then Xi := one of these nodes, mark G2 as tried for i
                if self.edgesleft > self.bestedgesleft or (runall and self.edgesleft >= self.bestedgesleft):
                    if self.G1node == self.g1_nodes.index(self.g1_nodes[-1]): # if the algorithm reached the end of the vertex list 
                        
                        #print "found MCES of size >> ",self.edgesleft,"in iter >> ",iter_cnt
                        mces_cnt+=1  
                        current = copy.deepcopy(self.current_mapping)
                        self.collect_best_assignments(current)
                        self.bestedgesleft = copy.deepcopy(self.edgesleft)    # dynamically assign bestedges score to edgesleft score every time MCS is found 

                        # routing for a look up of leftover assignment options for the final vertex removed from here [26/10/16] - Iva
                        # indicate that MCS has been found --> this will influence the treeorder appending
                        #self.bestedgesleft = copy.deepcopy(self.edgesleft)    # dynamically assign bestedges score to edgesleft score every time MCS is found 
                        timeittook = time.time() - self.starttime
                        #self.output.write('%s\n'%(timeittook))
                        self.output.write('%s\n'%(self.edgesleft))
                                               
                        for key,value in self.current_mapping.items():                          
                            self.output.write('%s\t%s\n'%(self.rG1indices[key],self.rG2indices[value]))
                        self.output.write('\n')
                        self.output.flush()
                        self.store_medges_edgesleft()
                        elapsed=time.time()-self.starttime                        
                        if (n_mces and mces_cnt >= n_mces) or (time_check and elapsed>maximum_time):
                            self.output.close()
                            final_assign_solutions = self.get_final_assignments()
                            return self.bestedgesleft,final_assign_solutions
                    else:  
                        self.store_medges_edgesleft()
                        self.G1node= self.G1node + 1
                        self.storage[self.G1node][2]=[] # mark all nodes as untried
                else:
                    # bring back unchanges medges for next iteration / G2 option checking
                    self.medges = medges_tmp
            else:       
                # if all matching options for current node have been tried
                # move a node up in a tree search              
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
                print "found MCES of size >> ",self.edgesleft,"in iter >> ",iter_cnt
                break
        self.output.close()
        final_assign_solutions = self.get_final_assignments()
        return self.bestedgesleft,final_assign_solutions
