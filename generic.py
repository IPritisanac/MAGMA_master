import os
import matplotlib as mpl
import matplotlib.pyplot as plt

class GenericMethods(object):
    #    methods of this class are generally useful and are therefore inherited by the other classes    
    #    returns overlap of the two lists/arrays and the reminder of either
    def Overlap(self,lista,listb):
        import collections
        a_multiset = collections.Counter(lista)
        b_multiset = collections.Counter(listb)
        overlap = list((a_multiset & b_multiset).elements())
        
        overlap = list((a_multiset & b_multiset).elements())
        a_remainder = list((a_multiset - b_multiset).elements())
        b_remainder = list((b_multiset - a_multiset).elements())
        return overlap,a_remainder,b_remainder
    
    def MakeLV_L(self,nodes,adjacency):
        new_node_labels = []
        for node in nodes:
            if node[-1] == "V":      
                new_node_labels.append(node[:-1] + "L")
            else:
                new_node_labels.append(node)
        
        new_adjacency = []
        for adj in adjacency:
            new_adj = []
            for a in adj:
                if a[-1] == "V":
                    new_adj.append(a[:-1] + "L")
                else:
                    new_adj.append(a)
            new_adjacency.append(new_adj)
        return new_node_labels, new_adjacency              
        
    def GetNodesAdjacency(self,contact_dict,lig_flag=False):
        if lig_flag:
            # collect node list and adjacency for all non ligand atoms
            node_list = [key for key in sorted(contact_dict.keys())]
            adjecancy = [list(set(contact_dict[node_list[n]])) for n in range(len(node_list))]
        else:
            node_list = [key for key in sorted(contact_dict.keys(),key = lambda x:int(x[:-1]))]
            adjecancy = [list(set(contact_dict[n])) for n in node_list]
        return node_list,adjecancy     
    
    def GetUnsortedNodesAdjacency(self,contact_dict):
        node_list = [key for key in contact_dict.keys()]
        adjecancy = [list(set(contact_dict[n])) for n in node_list]
        return node_list,adjecancy    

    def GetNodesAdjacencyDirect(self,contact_dict):
        node_list = [key for key in sorted(contact_dict.keys(),key = lambda x:int(x[1:]))]
        adjecancy = [list(set(contact_dict[n])) for n in node_list]
        return node_list,adjecancy
    
    def GetConnectionDictFromAdj(self,nodelist,adjacency):
        # given node list and adjacency
        # returns dictionary where keys are nodes and values are nodes connected to the key (node)
        conn_dict = {}
        for n in range(len(nodelist)):
            conn_dict.setdefault(nodelist[n],adjacency[n])
        return conn_dict       
    
    def SortDictType(self,dictionary):
        # given dictionary where keys are strings and the last character in the string indicated type
        # split the dictionary according to types
        types_dictionary = {}
        for key in dictionary.keys():
            types_dictionary.setdefault(key[-1],{})[key] = dictionary[key]
        return types_dictionary

    def RemoveFromPriority(self,remove,init_prior):
        new_init_prior = {key:value for key,value in init_prior.items() if key not in remove}
        return new_init_prior  
    
    def KnownAssignments(self,assignment_file):
        if os.path.isfile(assignment_file):
            known_assignments = {}
            fo = open(assignment_file,"r")
            for line in fo:
                stripped = line.strip()
                splitted = stripped.split()
                if len(splitted) == 2:
                    known_assignments.setdefault(splitted[0],[]).append(splitted[1])
                else:
                    continue
        else:
            print "No assignments known from before ..."
            known_assignments = {}
        return known_assignments
    
    def FixCandidatesDict(self,assign_candidates_dict,known_assignments):
        for key,value in assign_candidates_dict.items():
            if key in known_assignments.keys():
                assign_candidates_dict[key] = known_assignments[key]
        return assign_candidates_dict
        
    def JoinLigandConnections(self,protein_dict,ligand_dict):
        # join the values of the common keys between ligand and protein dictionaries
        # add the keys that do not exist in the protein dictionary to the protein dictionary        
        for key,value in ligand_dict.items():
            if key in protein_dict.keys():
                protein_dict[key]= list(set(protein_dict[key]+value))
            else:
                protein_dict[key]=value
        return protein_dict       
    
    def PlotAssignmentScoreMatrix(self,score_matrix,Npeaks,Natoms,label):
        
        cdict = {'red': ((0.0, 0.0, 0.0),
                        (0.5, 1.0, 0.0),
                        (1.0, 0.5, 1.0)),
                'green': ((0.0, 0.0, 0.0),
                          (0.5, 1.0, 0.0),
                          (1.0, 0.5, 1.0)),
                'blue': ((0.0, 0.0, 0.0),
                        (0.5, 1.0, 0.0),
                        (1.0, 0.5, 1.0))}
            
            
        cdict = {'red':   [(0.0,  1.0, 1.0),
                            (0.5,  1.0, 1.0),
                            (1.0,  0.0, 0.0)],
    
                    'green': [(0.0,  1.0, 1.0),
                              (0.25, 1.0, 1.0),
                              (0.75, 0.0, 0.0),
                              (1.0,  0.0, 0.0)],
    
                    'blue':  [(0.0,  1.0, 1.0),
                              (0.5,  0.0, 0.0),
                              (1.0,  0.0, 0.0)]}
        
        my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
        
        plt.pcolor(score_matrix, cmap=my_cmap)
        #plt.pcolor(self.pattern_score_matrix)
        plt.xticks(range(0,Npeaks,5))
        plt.yticks(range(0,Natoms,5))
        plt.title('Peak_residue mapping_'+str(label))
        plt.xlabel('Atom ID')
        plt.ylabel('Peak ID')
        plt.grid(True)
        plt.colorbar()
        plt.show()
        
 