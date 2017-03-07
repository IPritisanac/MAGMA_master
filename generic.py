"""
__author__: Iva Pritisanac

some generic methods used by several other classes (PDBData, NMRData, McGregor, MCES_PY,etc)
"""

import os
#import matplotlib as mpl
#import matplotlib.pyplot as plt

class GenericMethods(object):
    #    methods of this class are generally useful and are therefore inherited by the other classes    
    #    returns overlap of the two lists/arrays and the reminder of either
    def extract_overlap(self,lista,listb):
        import collections
        a_multiset = collections.Counter(lista)
        b_multiset = collections.Counter(listb)
        #overlap = list((a_multiset & b_multiset).elements())
        
        overlap = list((a_multiset & b_multiset).elements())
        a_remainder = list((a_multiset - b_multiset).elements())
        b_remainder = list((b_multiset - a_multiset).elements())
        return overlap,a_remainder,b_remainder

    def get_nodes_adjacency(self,contact_dict,lig_flag=False,merge_proRS=True): #  lig_flag and merge_proRS -- not needed
        
        node_list = [key for key in contact_dict.keys()] #,key = lambda x:int(x[:-1]))]
        adjecancy = [list(set(contact_dict[node_list[n]])) for n in range(len(node_list))]
        
        return node_list,adjecancy   
    
    def get_unsorted_nodes_adjacency(self,contact_dict):
        node_list = [key for key in contact_dict.keys()]
        adjecancy = [list(set(contact_dict[n])) for n in node_list]
        return node_list,adjecancy    

    
    def get_conns_dict_from_adj(self,nodelist,adjacency):
        # given node list and adjacency
        # returns dictionary where keys are nodes and values are nodes connected to the key (node)
        conn_dict = {}
        for n in range(len(nodelist)):
            conn_dict.setdefault(nodelist[n],adjacency[n])
        return conn_dict       
    
    def sort_dict_type(self,dictionary):
        # given dictionary where keys are strings and the last character in the string indicated type
        # split the dictionary according to types
        types_dictionary = {}
        for key in dictionary.keys():
            types_dictionary.setdefault(key[-1],{})[key] = dictionary[key]
        return types_dictionary
        
        
    def join_ligand_connections(self,protein_dict,ligand_dict):
        # join the values of the common keys between ligand and protein dictionaries
        # add the keys that do not exist in the protein dictionary to the protein dictionary        
        for key,value in ligand_dict.items():
            if key in protein_dict.keys():
                protein_dict[key]= list(set(protein_dict[key]+value))
            else:
                protein_dict[key]=value
        return protein_dict       