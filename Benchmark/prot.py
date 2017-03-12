#!/usr/bin/python

##########################################
#Get protparam details from a pdb file
#24th Oct 2016.
#A.Baldwin

import os,sys


from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from Bio.Data import IUPACData 


class entry():
    def __init__(self,test):
        self.pdb=test[0]
        self.dg=test[1]
        self.dgerr=test[2]
        self.nsb=test[3]
        self.nhb=test[4]
        self.oli=test[5]
        self.tds=test[6]
        self.newbury=test[7]
        self.bury=test[8]

def GetPDBs():
    inny=open('data.run')
    cnt=0
    pdb=[]
    dat=[]
    for line in inny.readlines():
        cnt+=1
        #if(cnt!=1):
        test=line.split()
        if(len(test)!=0):
                inst=entry(test)
                pdb.append(test[0])
                dat.append(inst)
    return pdb,dat

def AddIfNew(arr,test):
    for i in range(len(arr)):
        if(arr[i]==test):
            print 'found!'
            return
    arr.append(test)

#######################################


#pdbs,dat=GetPDBs()

inny=open('data.run')
for line in inny.readlines():
    pdb=line.split()[0]
    print 
    print 'analysing',pdb
    p=PDBParser(PERMISSIVE=1)
    structure_id=pdb
    #gunzip
    filename=pdb
    s=p.get_structure(structure_id,filename)

    ppb=PPBuilder()
    seqs=[]
    for pp in ppb.build_peptides(s):
        AddIfNew(seqs,str(pp.get_sequence()))

    anal=PA(seqs[0])
    per=anal.get_amino_acids_percent()
    por=anal.count_amino_acids()
    print 'MW:           ',anal.molecular_weight()
    print 'Aromaticity:  ',anal.aromaticity()
    print 'GRAVY:        ',anal.gravy()
    print 'InstabIndex:  ',anal.instability_index()
    print 'PI:           ',anal.isoelectric_point()
    print 'SSS:          ',anal.secondary_structure_fraction()
    print 'I:            ',por['I']
    print 'L:            ',por['L']
    print 'V:            ',por['V']
    print 'ILV methyls:  ',por['I']+por['V']*2+por['L']*2.
    print por
    print



#print pdbs,dat
