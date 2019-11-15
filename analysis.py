#!/usr/bin/python

"""
__author__: Andy Baldwin, Oct 2016
functions to analyse outputs of MAGMA calculations
"""
import os,sys

#get max edges score (last integer on a single line)
def get_score(infile):
    
    inny=open(infile)
    for line in inny.readlines():
        test=line.split()
        if(len(test)==1):
            try:
                val=int(test[0])    
            except:
                pass
                #print 'could not parse this',test[0]
    try:
        return val
    except:
        print('Could not find val in ',infile)
        sys.exit(100)


def add_file(result_dict,infile):
    
    test=infile.split('vf2')
    if(len(test)>1): #treat vf2 files differently
        lim='vf2'
        n_solutions=0
        results=[]
        fin = open(infile,"r")
        for line in fin:
            splitted = line.split()
            if(len(splitted)==2):
                key=splitted[0]
                val=splitted[1]
                if key not in result_dict.keys():
                    result_dict[key]=[]
                if val not in result_dict[key]:
                    result_dict[key].append(val)
    else:    #treat mces  .txt files differently
        lim=get_score(infile)
        print('The maximum NOE score is :',lim)
        tot_solutions=0
        n_solutions=0
        results=[]
        fin = open(infile,"r")
        for line in fin:
            splitted = line.split()
            if(len(splitted)==1):
                tot_solutions+=0.5 #two entries with line.split=1
                ass={}
                if(int(splitted[0])==lim):
                    n_solutions+=1
                    collect=True
                else:
                    collect=False
            elif(len(splitted)==0):
                if(len(ass.keys())!=0):
                    if(collect==True):
                        results.append(ass)
                    ass={}
                    collect=False
            elif(len(splitted)==2):
                ass.setdefault(splitted[0],[])
                if splitted[1] not in ass[splitted[0]]:
                    ass[splitted[0]].append(splitted[1])   
            #ass.append((splitted[0],splitted[1]))

        #print len(ass.keys()),collect
        if(len(ass.keys())!=0):
            if(collect==True):
                results.append(ass)

        fin.close()


        #print n_solutions,tot_solutions,len(results)
        for res in results:
            for key,vals in res.items():
                if key not in result_dict.keys():
                    result_dict[key]=[]
                for val in vals:
                    if val not in result_dict[key]:
                        result_dict[key].append(val)

    #return result_dict

def write_result_dict(result_dict,outfile):
    fout = open(outfile,"w")
    for k in sorted(result_dict,key=lambda k:len(result_dict[k])):
        #print k,
        fout.write("%s\t:"%(k))       
        for v in result_dict[k]:
            fout.write("%s\t"%(v))       
        fout.write('\n')
    fout.close()
    
def analyse_all(analdir):
    result_dict= {}
    files=os.listdir(analdir)
    for file in files:
        if(file[-4:]=='.txt'):
            test=file.split('new') #no 'new' files
            if(len(test)==1):
                tast=file.split('analysed') #no 'analysed' files
                if(len(tast)==1):
                    add_file(result_dict,analdir+os.sep+file)

    #incorr,corr=RunAnal(result_dict)                
    #print len(corr),corr,len(incorr),incorr
    return result_dict