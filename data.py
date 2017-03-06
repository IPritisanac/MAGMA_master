from generic import GenericMethods
import re,sys,os
from molecule import Molecule
from myparser import Parser
import numpy as np
import numpy.ma as ma
import scipy.spatial.distance as S
import matplotlib as mpl
import matplotlib.pyplot as plt

class ParseInput(Parser):
    #wrapper class that inherits all the methods of a class Parser
    #the class also inherits the constructor of the Parser
    ##Parser allows the extraction of keywords from an input text file 
    #input text file should contain a keyword per line, followed by one or more parameters (separated from the keyword by a tab)
    def __init__(self):
        
        super(ParseInput,self).__init__()
        self.add("use_pdb_file","pdb_flag","str","on")
        self.add('pdb_file','pdb_file','str',"")
        self.add('chains','chains','array str',[""])
        self.add('residues','residues','array str',[""])
        self.add('short_distance_threshold','short_distance_threshold','float',10.00)
        self.add('long_distance_threshold','long_distance_threshold','float',15.00)
        self.add('distance_type','distance_type','str','carbon')
        #self.add("rescore_threshold","rescore_thresh","float",8.0)
        self.add("distance_matrix","dist_matrix","str","")
        self.add("index_file","indx_file","str","")

        self.add('HMQC_list','HMQC_list','str',"")
        self.add('format_peak_list','format_peak_list','str',"")
        self.add('labelled_residues','labelled_residues','array str',[""])
        self.add('read_NOEs_from_file','read_NOEs_from_file','str',"")
        self.add('NOE_list','NOE_list','str',"")
        self.add('format_NOE_list','format_NOE_list','str',"")
        self.add('dimensions_NOEs','dimensions_NOEs','int',-1)
        self.add('read_extracted_NOEs','read_extracted_NOEs','str',"off")
        self.add('extracted_NOE_list','extracted_NOE_list','str',"")
        self.add('check_NOE_reciprocity','check_NOE_reciprocity','str',"on")
        self.add('check_NOE_intensity','check_NOE_intensity','str',"off")
        self.add('include_ligand_noes','include_ligand_noes','str','off')
        self.add('protein_ligand_noe_file','protein_ligand_noe_file','str',"")
        
        self.add("run_subgraph_isomorphism","check_subgraph_iso","str","on")
        self.add("priority_iterations","prior_iter","int",20)
        self.add("N_priority_MCES","n_prior_mces","int",3)
        self.add('output_file_McGregor','outfile_McGregor','str',"")
        self.add('output_file_vf2','outfile_vf2','str',"")               
        
        self.add('merge_NOEs','merge_NOEs','str',"")
        self.add('merge_proRS','merge_proRS','str',"")
        self.add('merge_LV_label','merge_LV_label','str',"")
        self.add("maximum_run_time", "max_run_time","int","")
        self.add('min_size','min_size','int',3)
        self.add('mode','mode','str',"complete_subgraphs")
        self.add('run_at','run_at','str',"short")
        self.add("mces_mode","mces_mode","str","all")

        #self.add("craicMode","craicMode","str","all")
        #self.add("craicFile","craicFile","str","all")
        self.add("strip_mode","strip_mode","str","on")

        self.add("include_ligand","include_ligand",'str',"off")
        self.add("ligand_chain","ligand_chain",'str',"")
        self.add("ligand_name","ligand_name",'str',"")
        self.add('ligand_atoms','ligand_atoms','array str',[""])
        self.add("ligand_distance_threshold","ligand_distance_threshold","int",5)

    def check_variables(self):
        
        if self.pdb_flag == "":
            raise Exception("in input file decide if pdb file is used or not!\n In the line use_pdb_file put 'off' if the file should not be used or 'on' if it should be used")
        elif self.pdb_flag == "off":
            self.pdb_flag = False
        elif self.pdb_flag == "on":
            self.pdb_flag = True
        else:
            raise Exception("For 'use_pdb_file' expected 'on' or 'off', got %s instead"%(self.pdb_flag))
        
        if self.pdb_file == "":
            raise Exception("pdb file must be defined!")

        if os.path.isfile(self.pdb_file) == False:
            raise Exception("pdb file does not exist or the path to the file is incorrect!")

        if self.dist_matrix and not os.path.isfile(self.dist_matrix):
            raise Exception("distance matrix is specified, but no file is provided!\nEither exclude distance_matrix line from an input file or define a file!")

        if self.indx_file and not os.path.isfile(self.indx_file):
            raise Exception("index file is specified, but no file is provided!\nEither exclude the line from an input file or define a file!")

        if self.dist_matrix and not self.indx_file:
            raise Exception("distance matrix is provided, but the index file of the matrix is not!")
        
        if len(self.chains)==0:          
            """
            //Default could be set instead of raising an exception
            """  
            raise Exception("no chains of the pdb file selected")
        
        if len(self.residues) == 0:
            """
            //Default could be set instead of raising an exception
            """
            raise Exception("no methyl residues of the pdb file selected")
        
        if self.short_distance_threshold == "":
            """
            //Default could be set instead of raising an exception
            """
            raise Exception("short distance threshold for the calculation is not provided ")
            
        if self.long_distance_threshold == "":
            """
            //Default could be set instead of raising an exception
            """
            raise Exception("long distance threshold for the calculation is not provided ")
        
        if (self.distance_type!="carbon") and (self.distance_type!="proton"):
            raise Exception("Distance type can be either 'carbon' or 'proton'. Instead got %s"%self.distance_type)
        
        if self.distance_type=="":
            raise Exception("Please define between which methyl atoms distance matrix should be computed: 'proton' or 'carbon'")
        
        #if self.rescore_thresh == "":
            #print "No rescore against different distance threshold chosen!"
            #raise Exception("No rescore against different distance threshold chosen!")
        
        if self.short_distance_threshold > self.long_distance_threshold:
            raise Exception("short and long distance thresholds might have been swapped!\ncheck short/long distance thresholds!")
        
        if self.HMQC_list == "":
            raise Exception("peak/HMQC list must be given!")
        
        if os.path.isfile(self.HMQC_list) == False:
            raise Exception("peak list file (HMQC peak list) does not exist or the path to the file is incorrect!")
        
        if len(self.labelled_residues) == 0:
            raise Exception("no methyl labelled residues selected!\n list residues with methyl labels!")
        
        if sorted(self.labelled_residues) != sorted(self.residues):
            raise Exception("methyl labels from PDB file do not match with the NMR input?\n make sure that list 'residues' and list 'labelled residues' match!")
        
        if self.format_peak_list not in ["sparky"]:
            raise Exception("only sparky input format of the HMQC peak list is currently supported!")
        
        if self.read_NOEs_from_file!= "on" and self.read_NOEs_from_file!="off":
            raise Exception("incorrect entry for 'read_NOE_from_file' ! Expected on or off; instead got %s"%self.read_NOEs_from_file)
                
        if self.read_NOEs_from_file == "on" and self.read_extracted_NOEs == "on":
            raise Exception("incorrect entries for 'read_NOE_from_file' and/or 'self.read_extracted_NOEs'! Only one of the two NOE modes can be on!")

        if self.read_NOEs_from_file == "on" and os.path.isfile(self.NOE_list)== False:
            raise Exception("'read_NOEs_from_file' is on, however 'NOE_list' file does not exist or the path to the file is incorrect!")
                    
        if self.read_NOEs_from_file == "on" and self.NOE_list == "":
            raise Exception("if read_NOEs_from_file is on, 'NOE_list' file must be given!")
        
        if os.path.isfile(self.NOE_list)== False:
            raise Exception("NOE file does not exist or the path to the file is incorrect!")
        
        if self.read_NOEs_from_file == "on":
            self.read_NOEs_from_file = True
        else:
            self.read_NOEs_from_file = False
        
        if self.read_NOEs_from_file and self.format_NOE_list not in ["sparky"]:
            raise Exception("only sparky input format of the NOE list is supported!")
        
        if self.dimensions_NOEs == -1:
            raise Exception("number of dimensions of the NOESY experiment not given!")
        
        if self.dimensions_NOEs != 3 and self.dimensions_NOEs != 4:
            raise Exception("incorrect number of dimensions in the NOE file given! Expected 3 or 4; instead got %s"%self.dimensions_NOEs)
        
        if self.read_extracted_NOEs != "off" and self.read_extracted_NOEs !="on":
            raise Exception("incorrect entry for 'read_NOE_directly' ! Expected on or off; instead got %s"%self.read_extracted_NOEs)
                    
        if self.read_extracted_NOEs == "on" and self.extracted_NOE_list == "":
            raise Exception("if 'read_extracted_NOEs' is on, 'extracted_NOE_list' file must be given!")
                       
        if self.read_extracted_NOEs == "on" and os.path.isfile(self.extracted_NOE_list) == False:
            raise Exception("'read_extracted_NOEs' is on, however 'extracted_NOE_list' file does not exist or the path to the file is incorrect!")
                       
        if self.read_extracted_NOEs == "on":
            self.read_extracted_NOEs = True
        else:
            self.read_extracted_NOEs = False
        
        if not self.read_extracted_NOEs and not self.read_NOEs_from_file:
            raise Exception("One of 'read_extracted_NOEs' or 'read_NOEs_from_file' needs to be on!")
                        
        if self.check_NOE_reciprocity != "" and self.check_NOE_reciprocity != "on" and self.check_NOE_reciprocity != "off":
            raise Exception("incorrect entry for 'check_NOE_reciprocity'; expected 'on' or 'off'; instead got %s"%self.check_NOE_reciprocity)
        
        if self.check_NOE_intensity != "" and self.check_NOE_intensity != "on" and self.check_NOE_intensity != "off":
            raise Exception("incorrect entry for 'check_NOE_intensity'; expected 'on' or 'off'; instead got %s"%self.check_NOE_intensity)
        
        if self.check_NOE_reciprocity == "on":
            self.check_NOE_reciprocity = True
        elif self.check_NOE_reciprocity == "off":
            self.check_NOE_reciprocity = False   
        
        if self.check_NOE_intensity == "on":
            self.check_NOE_intensity = True
        elif self.check_NOE_intensity == "off":
            """Default currently >> off"""
            self.check_NOE_intensity = False
        else:
            raise Exception("Incorrect entry for check_NOE_intensity, expected 'on' or 'off'. Instead got %"%self.check_NOE_intensity)
        
        if self.include_ligand_noes=="on":
            self.include_ligand_noes = True
        elif self.include_ligand_noes=="off":
            self.include_ligand_noes = False
        else:
            raise Exception("incorrect entry for 'include_ligand_noes'; expected 'on' or 'off'; instead got %s"%self.include_ligand_noes)
        
        if self.include_ligand_noes and not os.path.isfile(self.protein_ligand_noe_file):
            raise Exception("'include_ligand_noes' is on, however 'protein_ligand_noe_file' file does not exist or the path to the file is incorrect!")
                    
        if self.outfile_McGregor == "":
            self.outfile_McGregor = "output_MCES_run.txt"
            #raise Exception("file for results of McGregor MCES search does not exist or the path to the file is incorrect!")
        
        if self.outfile_vf2 == "":
            self.outfile_vf2 = "output_subgraph_isomorphism_vf2_run.txt"
            #raise Exception("file for results of vf2 subgraph isomorphism search does not exist or the path to the file is incorrect!")
       
        if self.check_subgraph_iso=="on":
            self.check_subgraph_iso = True
        elif self.check_subgraph_iso=="off":
            self.check_subgraph_iso = False
        else:
            raise Exception("incorrect entry for 'check_subgraph_iso'; expected 'on' or 'off'; instead got %s"%self.check_subgraph_iso)  
                    
        if self.merge_NOEs == "on":
            self.merge_NOEs = True
        elif self.merge_NOEs == "off":
            self.merge_NOEs = False
        else:
            raise Exception("incorrect entry for 'merge_NOEs'; expected 'on' or 'off'; instead got %s"%self.merge_NOEs)
        
        if self.merge_proRS == "on":
            self.merge_proRS = True
        elif self.merge_proRS == "off":
            self.merge_proRS = False
        else:
            raise Exception("incorrect entry for 'merge_proR_S'; expected 'on' or 'off'; instead got %s"%self.merge_proR_S)
                        
        if self.merge_LV_label == "on":
            self.merge_LV_label = True
        elif self.merge_LV_label == "off":
            self.merge_LV_label = False
        else:
            raise Exception("incorrect entry for 'merge_LV_label'; expected 'on' or 'off', instead got %s"%(self.merge_LV_label))
                                
        if  not isinstance(self.min_size,int):
            raise Exception("incorrect entry for 'min_size'; integer expected, instead got %s"%(self.min_size))
                              
        if self.mode not in ["all_subgraphs", "connected_subgraphs","residue_subgraphs"]:
            raise Exception("mode keyword should be all_subgraphs, connected_subgraphs or residue_subgraphs. Instead found: %s"%self.mode)

        """
        //SET proper defaults and checks for this//        
        if self.strip_mode == "on":   # is set to on or a default applies
            self.strip_mode = "on"
        elif self.strip_mode == "off":
            self.strip_mode = None
        else:
            raise Exception("incorrect entry for 'strip_mode'; expected 'on' or 'off', instead got %s"%(self.strip_mode))
        """
        if self.run_at == "":
            raise Exception("distance threshold for the calculation must be defined!/n set 'run_at' to either 'short' or 'long'")
        
        if (self.run_at!="short") and (self.run_at!="long"):
            raise Exception("incorrect entry for'run_at'; expected 'short' or 'long', instead got %s"%(self.run_at))
        
        if self.prior_iter == "":
            raise Exception("number of iterations for priority calculation must be defined!")
        
        if self.n_prior_mces == "":
            raise Exception("set a number for MCES to be generated in priority setting step!")          
        
        if self.max_run_time == "":
            raise Exception("set an integer value or None for the maximum_run_time")
                
        if self.mces_mode == "":
            raise Exception("running mode for MCES algorithm must be defined!/n set 'mces_mode' to either 'all' or 'one'")
           
        if self.include_ligand == "on":
            self.include_ligand = True
        elif self.include_ligand == "off":
            self.include_ligand = False
        else:
            raise Exception("incorrect entry for 'include_ligand'; expected 'on' or 'off', instead got %s"%(self.include_ligand))
        
        if self.include_ligand and (self.ligand_chain==""):
            raise Exception("include_ligand is 'on', but ligand_chain is not given! Please enter chain of the ligand; check the pdb file")

        if self.include_ligand and (self.ligand_name==""):
            raise Exception("include_ligand is 'on', but ligand_name is not defined! Please enter name/ID of the ligand; check the pdb file")

        if self.include_ligand and (len(self.ligand_atoms)==0):
            raise Exception("include_ligand is 'on', but ligand_atoms are not given! Please enter atoms of the ligand of interested; separated by a space")
        
        if self.include_ligand and (self.ligand_distance_threshold == ""):
            raise Exception("include_ligand is 'on', but ligand_distance_threshold is not given! Please enter integer value for ligand_distance_threshold")
                                                              
class PDBData(GenericMethods):
    
    def __init__(self,input_file,input_chains,input_residues,dist_between,short_thresh,long_thresh,merge_proRS=True,ligand=False,lig_chain=None,lig_name=None,lig_atoms=None):
        
        self.input_file = input_file
        self.input_chains = input_chains
        self.input_residues = input_residues
        self.dist_between = dist_between
        self.short_thresh = short_thresh
        self.long_thresh = long_thresh
        self.merge_proRS = merge_proRS
        
        self.ligand = ligand
        self.lig_chain = lig_chain
        self.lig_name = lig_name
        self.lig_atoms = lig_atoms   

        self.methyl_proton_atom_names = {"ILE":["HD11","HD12","HD13"],"LEU":["HD11","HD12","HD13","HD21","HD22","HD23"],"VAL":["HG11","HG12","HG13","HG21","HG22","HG23"],"THR":["HG21","HG22","HG23"],"ALA":["HB1","HB2","HB3"],"MET":["HE1","HE2","HE3"]}        
        # start empty data structures
        self.all_res_information = {}
        
    def CheckInputFiles(self):
        #   method checks that all required input files exist    
        pass
    
    def CheckInputFormats(self):
        #   method checks that all input files are in correct format 
        pass
           
    def PDBHandle(self):
        #    check whether pdb file exists and create object of class Molecule
        #    return this object or exit the program        
        self.input_chains= list(self.input_chains)        
        M=Molecule()   
        try:
            M.import_pdb(self.input_file)
        except:
            raise Exception('ERROR: file %s not found!'% self.input_file)
            #print "File %s not found!"  % self.input_file
            #print 'Exiting program ...'
            #sys.exit()
        return M       
        
    def ExtractPositionsIndices(self,flag):  
        # extract positions and indices of labelled atoms from a pdb file        
        # in case we are not reading the information in from the pdb file - exit here, and read in the information from an external file in a different format

        if not flag:    # if information is read in from a different, suitably formated file
            return None,None
        
        mol = self.PDBHandle()        
        print "> protein loaded, loading aminoacids!"
        
        if self.dist_between == "proton":
            all_pos = []
            all_info = []
               
        if 'ILE' in self.input_residues:
            #get position, indices and information for isoleucine atoms
            if self.dist_between == "carbon":
                pos_ile,idx_ile=mol.atomselect(self.input_chains,"ILE","CD1",use_resname=True, get_index=True)
                info_ile=mol.properties['data'][idx_ile,2:6]    
                print ">> isoleucine done!",pos_ile.shape
                self.all_res_information.setdefault('ILE',(pos_ile,info_ile))
                
            if self.dist_between == "proton":
                for patom in self.methyl_proton_atom_names["ILE"]:
                    pos,idx = mol.atomselect(self.input_chains,"ILE",str(patom),use_resname=True,get_index=True)
                    info = mol.properties['data'][idx,2:6]
                    all_pos.extend(pos)
                    all_info.extend(info)
                print ">> isoleucine done!", len(all_pos)
            
        if 'ALA' in self.input_residues:
        #get position, indices and information for alanine atoms
            if self.dist_between == "carbon":
                pos_ala,idx_ala=mol.atomselect(self.input_chains,"ALA","CB",use_resname=True, get_index=True)
                info_ala=mol.properties['data'][idx_ala,2:6]
                print ">> alanine done!",pos_ala.shape
                self.all_res_information.setdefault('ALA',(pos_ala,info_ala))
            
            if self.dist_between == "proton":
                for patom in self.methyl_proton_atom_names["ALA"]:
                    pos,idx = mol.atomselect(self.input_chains,"ALA",str(patom),use_resname=True,get_index=True)
                    info = mol.properties['data'][idx,2:6]
                    all_pos.extend(pos)
                    all_info.extend(info)
                print ">> alanine done!", len(all_pos)
        
        if 'THR' in self.input_residues:                
            if self.dist_between == "carbon":          
                #get position, indices and information for threonine atoms
                pos_thr,idx_thr=mol.atomselect(self.input_chains,"THR","CG2",use_resname=True, get_index=True)
                info_thr=mol.properties['data'][idx_thr,2:6]
                print ">> threonine done!"
                self.all_res_information.setdefault('THR',(pos_thr,info_thr))
                
            if self.dist_between == "proton":
                for patom in self.methyl_proton_atom_names["THR"]:
                    pos,idx = mol.atomselect(self.input_chains,"THR",str(patom),use_resname=True,get_index=True)
                    info = mol.properties['data'][idx,2:6]
                    all_pos.extend(pos)
                    all_info.extend(info)
                print ">> threonine done!", len(all_pos)
                        
        if 'MET' in self.input_residues:               
            #get position, indices and information for methionine atoms
            if self.dist_between == "carbon": 
                pos_met,idx_met=mol.atomselect(self.input_chains,"MET","CE",use_resname=True, get_index=True)
                info_met=mol.properties['data'][idx_met,2:6]
                print ">> methionine done!"
                self.all_res_information.setdefault('MET',(pos_met,info_met))
                
            if self.dist_between == "proton":
                for patom in self.methyl_proton_atom_names["MET"]:
                    pos,idx = mol.atomselect(self.input_chains,"MET",str(patom),use_resname=True,get_index=True)
                    info = mol.properties['data'][idx,2:6]
                    all_pos.extend(pos)
                    all_info.extend(info)
                print ">> methionine done!", len(all_pos)
        
        if 'VAL' in self.input_residues: 
            #get average position, indices and information for valine atoms
            if self.dist_between == "carbon": 
                if self.merge_proRS:
                    idx_val=mol.atomselect(self.input_chains,"VAL","CA",use_resname=True, get_index=True)[1]
                    info_val=mol.properties['data'][idx_val,2:6]
                    print ">> merging ProR and ProS carbons to pseudoatom"
                    pos_val=self.GetAveragePos(mol,idx_val,["CG1","CG2"])
                    print ">> valine done!",pos_val.shape
                    self.all_res_information.setdefault('VAL',(pos_val,info_val))
                else:
                    print "Treating each methyl group of Val separately"
                    pos_val1,idx_val1=mol.atomselect(self.input_chains,"VAL","CG1",use_resname=True, get_index=True)
                    pos_val2,idx_val2=mol.atomselect(self.input_chains,"VAL","CG2",use_resname=True, get_index=True)
                    info_val1=mol.properties['data'][idx_val1,2:6]
                    info_val2=mol.properties['data'][idx_val2,2:6]
                    print ">> %s valine methyls found!"%np.concatenate((pos_val1,pos_val2)).shape[0] 
                    print ">> valine done!"
                    self.all_res_information.setdefault("VAL",(np.concatenate((pos_val1,pos_val2)),np.concatenate((info_val1,info_val2))))
                
            if self.dist_between == "proton":
                if self.merge_proRS:
                    for patom in self.methyl_proton_atom_names["VAL"]:
                        pos,idx = mol.atomselect(self.input_chains,"VAL",str(patom),use_resname=True,get_index=True)
                        info = mol.properties['data'][idx,2:6]
                        all_pos.extend(pos)
                        all_info.extend(info)
                    print ">> valine done!", len(all_pos)
                else:
                    print "Currently unsupported combination of parameters: distance_type 'proton' and merge_proRS 'off'; Change one of the parameters"
                    sys.exit(0)
            
        if 'LEU' in self.input_residues:        
            #get average position, indices and information for leucine atoms
            if self.dist_between == "carbon":
                if self.merge_proRS:
                    print "Treating each methyl group of Leu separately"           
                    idx_leu=mol.atomselect(self.input_chains,"LEU","CA",use_resname=True, get_index=True)[1]
                    info_leu=mol.properties['data'][idx_leu,2:6]
                    print ">> merging ProR and ProS carbons to pseudoatom"
                    pos_leu=self.GetAveragePos(mol,idx_leu,["CD1","CD2"])
                    print ">> leucine done!",pos_leu.shape
                    self.all_res_information.setdefault('LEU',(pos_leu,info_leu))
                else:
                    pos_leu1,idx_leu1=mol.atomselect(self.input_chains,"LEU","CD1",use_resname=True, get_index=True)
                    pos_leu2,idx_leu2=mol.atomselect(self.input_chains,"LEU","CD2",use_resname=True, get_index=True)
                    info_leu1=mol.properties['data'][idx_leu1,2:6]
                    info_leu2=mol.properties['data'][idx_leu2,2:6]
                    print ">> %s leucine methyls found!"%np.concatenate((pos_leu1,pos_leu2)).shape[0]
                    print ">> leucine done!"        
                    self.all_res_information.setdefault("LEU",(np.concatenate((pos_leu1,pos_leu2)),np.concatenate((info_leu1,info_leu2))))

            if self.dist_between == "proton":
                if self.merge_proRS:
                    for patom in self.methyl_proton_atom_names["LEU"]:
                        pos,idx = mol.atomselect(self.input_chains,"LEU",str(patom),use_resname=True,get_index=True)
                        info = mol.properties['data'][idx,2:6]
                        all_pos.extend(pos)
                        all_info.extend(info)
                    print ">> leucine done!", len(all_pos)
                else:
                    print "Currently unsupported combination of parameters: distance_type 'proton' and merge_proRS 'off'; Change one of the parameters"
                    sys.exit(0)
                                    
        if (self.dist_between == "carbon") and (len(self.all_res_information.keys()) == 0):
            print "No residues selected!"
            print "Check input of residues (>> PDB part of the input file)"
            
        if (self.dist_between == "carbon") and (len(self.all_res_information.keys()) == 1):
            print "Only 1 residue type methyl labelled >> ", self.all_res_information.keys()
            residue = self.all_res_information.keys()[0]
            all_pos = self.all_res_information[residue][0]
            all_info = self.all_res_information[residue][1]
            
        if (self.dist_between == "carbon") and (len(self.all_res_information.keys()) > 1):
            # initialize numpy arrays of right dimensions
            first_residue = self.all_res_information.keys()[0] 
            all_pos = self.all_res_information[first_residue][0]
            all_info = self.all_res_information[first_residue][1]
            for r in range(len(self.all_res_information.keys())-1):
                # append the rest
                res = self.all_res_information.keys()[r+1]
                all_pos=np.concatenate((all_pos,self.all_res_information[res][0]))
                all_info=np.concatenate((all_info,self.all_res_information[res][1]))
        
        if self.ligand:
            #    if atoms from ligand are connected to the protein residues through an NOE;
            #    and these NOEs will be included to the structure graph as edges
            print "\n>> adding edges from the ligand to the protein\n"
            for lig_atom in self.lig_atoms:
                pos_lig,idx_lig=mol.atomselect(str(self.lig_chain),str(self.lig_name),str(lig_atom),use_resname=True,get_index=True)
                info_lig=mol.properties['data'][idx_lig,2:6] 
                all_pos = np.concatenate((all_pos,pos_lig))
                all_info = np.concatenate((all_info,info_lig))
        
        all_pos = np.array(all_pos)
        all_info = np.array(all_info)
        
        return all_pos,all_info,all_info[:,1:]
    
    def GetAveragePos(self,mol,idx,atoms):
        #compute average position for atoms part of a residue identified by a specific index
        # M=molecule, idx=list of indices of atoms (in different resid), atoms=2 atoms to look for in residues of interest)    
        pos_tmp=[]
        for i in idx:
            pts,res_idxs=mol.same_residue(i,get_index=True)
            idx_subset=np.logical_or(mol.properties['data'][res_idxs,2]==atoms[0],mol.properties['data'][res_idxs,2]==atoms[1])
            pos_tmp.append(np.mean(pts[idx_subset],axis=0))
        return np.array(pos_tmp)
                
    def LigandInteractionNetwork(self,positions,info,lig_thresh):
        # if connections to ligand atoms are a part of the structure graph
        # do additional filtering of connections, based on the ligand specific threshold
        dist_mat=S.cdist(positions,positions)
        masked_distances=ma.MaskedArray(dist_mat, mask=np.logical_or(dist_mat>lig_thresh,dist_mat==0))
        #find couples in contact
        couples=np.array(np.where(np.logical_and(dist_mat<lig_thresh,dist_mat>0))).T                
        """
        check if ligand atom
        extract couples with ligand atom indicated
        """
        self.lig_atoms = list(self.lig_atoms)
        
        t=[]
        for c in couples:
            if c[0]>c[1]:
                continue
                        
            if (info[c[0]][1]==str(self.lig_name)) and (info[c[1]][1]==str(self.lig_name)):
                # if both ends of an edge are ligand atoms
                # skip the couple
                continue
            
            if (info[c[0]][1]==str(self.lig_name) and info[c[0]][0] in self.lig_atoms) or (info[c[1]][1]==str(self.lig_name) and info[c[1]][0] in self.lig_atoms):
                # if one of the edge ends it a ligand atom
                if info[c[0]][1]==str(self.lig_name):
                    name1 = info[c[0]][1]+info[c[0]][0][:-1]
                else:
                    name1=''.join(info[c[0]][1:4:2])
                    
                if info[c[1]][1]==str(self.lig_name):
                    name2 = info[c[1]][1]+info[c[1]][0][:-1]
                else:
                    name2=''.join(info[c[1]][1:4:2])
                
                #print info[c[0]], ">>", info[c[1]]
                #print name1,">>", name2
            
                if name1==name2:
                    continue            
                names=np.array([name1,name2])
                order=np.argsort(names)
                label='_'.join(names[order])
                t.append([label, c[order[0]],c[order[1]], dist_mat[c[order[0]],c[order[1]]]])
        
        test_table=np.array(t)        
        #find couple at shortest distance, if more than one occurrence
        newdata=[]
        for i in np.unique(test_table[:,0]):
            lines=test_table[test_table[:,0]==i]
            bestpos=np.argmin(lines[:,3].astype(float))
            newdata.append(lines[bestpos])
        newtable=np.array(newdata)
        return newtable,info,dist_mat,masked_distances
    
    def InteractionNetwork(self,positions,info,protein_thresh,posfile=None,indxfile=None):
        #    compute a distance map
        #    ignore doubles (double edges)
        #    ignore self connection (loops)
        #    find connection at shortest distance if more than one occurrence of the link between 2 atoms (in e.g. homodimers or higher order oligomers)
        #    in case distance matrix is being imported posfile and indxfile will be provided
        #    read those in ... 
        
        #red_info = info[:,1:]
        
        if posfile and indxfile:
            dist_mat = np.loadtxt(posfile)
            #dist_mat=S.cdist(positions,positions)
            info = np.loadtxt(indxfile,dtype=str)
        else:
            dist_mat=S.cdist(positions,positions)
            info = info           
        
        masked_distances=ma.MaskedArray(dist_mat, mask=np.logical_or(dist_mat>protein_thresh,dist_mat==0))
        #find couples in contact
        couples=np.array(np.where(np.logical_and(dist_mat<protein_thresh,dist_mat>0))).T
                
        t=[]
        for c in couples:
            if c[0]>c[1]:
                continue
                        
            if self.ligand and (info[c[0]][0]==str(self.lig_name) or info[c[1]][0]==str(self.lig_name)):
                # skip all ligand connections
                continue
            
            name1=''.join(info[c[0]][0:3:2])
            name2=''.join(info[c[1]][0:3:2])
                        
            if name1==name2:
                continue
        
            names=np.array([name1,name2])
            order=np.argsort(names)
            label='_'.join(names[order])
            t.append([label, c[order[0]],c[order[1]], dist_mat[c[order[0]],c[order[1]]]])

        test_table=np.array(t)        
        #find couple at shortest distance, if more than one occurrence
        newdata=[]
        for i in np.unique(test_table[:,0]):
            lines=test_table[test_table[:,0]==i]
            bestpos=np.argmin(lines[:,3].astype(float))
            newdata.append(lines[bestpos])
        newtable=np.array(newdata)
        
        return newtable,info,dist_mat,masked_distances
    
    def WriteIndices(self,atom_info,filename):
        fout = open(filename,"w")
        for i in range(len(atom_info)):
            #print i,atom_info[i][0]
            fout.write("%s\t%s\t%s\t%s\n"%(i,atom_info[i][0],atom_info[i][1],atom_info[i][2]))
        fout.close()
    
    def GetConnectionDict(self,test_table,all_atoms,ligand_conns=False):
        # given the table of pairs (atoms connected through space within distance threshold >> see self.InteractionNetwork
        # extract the pairs into dictionary where key is an atom and value is a list of atoms key is connected to
        # in case some atoms from all_atoms are not paired with any other one
        # add those as keys with empty values in the dictionary
        # return connection dictionary
        conn_dict_pdb = {}        
        for row in test_table:
            conn = row[0].split('_')
                        
            if ligand_conns and (conn[0][:3]==str(self.lig_name)):
                conn1 = conn[0][3:]+"X"
            else:
                conn1 = conn[0][3:] + conn[0][0]
            if ligand_conns and (conn[1][:3]==str(self.lig_name)):
                conn2 = conn[1][3:]+"X"
            else:
                conn2 = conn[1][3:] + conn[1][0]
            
            #print conn1,">>", conn2
            conn_dict_pdb.setdefault(conn1,[]).append(conn2)
            conn_dict_pdb.setdefault(conn2,[]).append(conn1)
       
        for key,value in conn_dict_pdb.items():
            newval = set(value)
            conn_dict_pdb[key] = list(newval)
        
        if self.ligand and ligand_conns:
            # check if all atoms are included in the dictionary, if not add those missing 
            ligatoms = [row[0][:-1]+"X" for row in all_atoms if str(row[1])==str(self.lig_name)]
            leftover_lig = list(set(ligatoms) - set(conn_dict_pdb.keys()))
            #print "ligand atom leftovers..."
            #print leftover_lig
            for g in leftover_lig:
                conn_dict_pdb[g] = []
        
        if self.ligand and not ligand_conns:    
            # check if all atoms are included in the dictionary, if not add those missing 
            all_atoms =  [row[2]+row[0][0] for row in all_atoms if str(row[0])!=self.lig_name]   
            leftover_atoms = list(set(all_atoms) - set(conn_dict_pdb.keys()))
            # add the missing ones to the dictionary with empty values (no connections)
            #print leftover_atoms
            for a in leftover_atoms:
                conn_dict_pdb[a] = []
                
        if not self.ligand and not ligand_conns:    
            all_atoms =  [row[2]+row[0][0] for row in all_atoms] 
            leftover_atoms = list(set(all_atoms) - set(conn_dict_pdb.keys()))
            # add the missing ones to the dictionary with empty values (no connections)
            for a in leftover_atoms:
                conn_dict_pdb[a] = []
                
        return conn_dict_pdb
    
    def GetConnectionDictDist(self,test_table,all_atoms):
        # given the table of pairs (atoms connected through space within distance threshold >> see self.InteractionNetwork
        # extract the pairs into dictionary where key is an atom and value is a list of atoms key is connected to
        # return connection dictionary in which for each connection, distance is specified as well
        conn_dict_pdb = {}
        for row in test_table:
            conn = row[0].split('_')
            conn1 = conn[0][3:] + conn[0][0]
            conn2 = conn[1][3:] + conn[1][0]
            conn_dict_pdb.setdefault(conn1,[]).append((conn2,float(row[3])))
            conn_dict_pdb.setdefault(conn2,[]).append((conn1,float(row[3])))
        for key,value in conn_dict_pdb.items():
            newval = set(value)
            conn_dict_pdb[key] = list(newval)
        return conn_dict_pdb    
    
    def WriteDistanceFile(self,conn_dict,filename):
        fout = open(filename,"w")
        for atom,conns in conn_dict.items():
            for c in conns:
                fout.write("%s\t%s\t%s\n"%(atom,c[0],c[1]))
        fout.close() 
                    
    def PrettyPlotContacts(self,table,positions,distance_matrix,mask_dist,thresh):
        #data for pretty plot here!!!
        d=[]
        for line in table:
          
            point1=positions[int(line[1])]
            point2=positions[int(line[2])]
            d.append([point1, point2, float(line[3])])
        #print d        
        #plot a pretty matrix
        print "> preparing matrix plot"
        length=np.arange(0,len(distance_matrix),1)
        [X,Y] =np.meshgrid(length,length)
        
        #fig=plt.figure(0)
        #ax=fig.add_subplot(1,1,1)
        #ax.pcolor(X,Y,mask_dist,cmap = mpl.cm.get_cmap("hot"),vmin=0, vmax= thresh)
        plt.figure()
        #ax=fig.add_subplot(1,1,1)
        plt.pcolor(X,Y,mask_dist,cmap = mpl.cm.get_cmap("hot"),vmin=0, vmax= thresh)
        plt.colorbar()
        plt.show()
            
    def TurnValToLeu(self,conn_dict):
        #    if valine and leucines are not discriminated
        #    turns all valine labels to leucine labels
        relabelled_conn_dict = {}
        for key,value in conn_dict.items():
            new_value = []
            if key[-1] == "V":
                new_key = key[:-1]+"L"
            else:
                new_key = key
            for v in value:
                if v[-1]=="V":
                    new_value.append(v[:-1]+"L")
                else:
                    new_value.append(v)
            relabelled_conn_dict.setdefault(new_key,new_value)
                
        return relabelled_conn_dict
        
    
class NMRData(GenericMethods):
    #the class extracts a set of vertices and edges from NMR data
    #the data needs to be provided in an appropriate format 
    #appropriate format for inter-methyl NOE data is two column file, separated by a tab, with entry single letter amino acid identifier+resonance number (e.g. I1    V2)
    #hmqc file is a three column file with resonance id in first column, 13C and 1H frequencies in 2nd and 3rd, respectively
    #hmqc file can be an empty file, as the methyl resonances with no NOEs are excluded from the data graph definition
    #different file formats will be supported in the futures, as indicated below    
    def __init__(self,input_peak_list,input_list_format,input_labelled_residues,input_noe,input_noe_list,input_noe_format,dimensions,read_extracted_NOEs,extracted_NOEs,ligand_noe_file,reciprocity=True,intensity=False,mergeLV=True,mergeLVlabel=False,ligand_noes=False):
        
        self.input_peak_list = input_peak_list
        self.input_list_format = input_list_format
        self.input_peak_residues = input_labelled_residues
        self.input_noe = input_noe
        self.input_noe_list = input_noe_list
        self.input_noe_format = input_noe_format
        self.noesy_dimension = dimensions
        
        self.read_extracted_NOEs = read_extracted_NOEs
        self.extracted_NOEs = extracted_NOEs
        self.check_recip = reciprocity
        self.check_intensity = intensity    
        self.merge = mergeLV
        self.merge_label_LV = mergeLVlabel
        
        self.ligand_noes = ligand_noes
        self.ligand_noe_file = ligand_noe_file
        
        self.methyl_IDs = [type[0] for type in self.input_peak_residues]
        self.peaks = []   
         
    def CheckInputFormats(self):
        """
        //TO DO//
        """
        #   method checks that all input files are in correct format 
        #    sets class variables that will be used by NOE reading and checking methods
        if self.input_noe_format == 'sparky':
            '''
            // TO DO    //
            >>    check dimensions
            >>    order of dimensions
            >>    consistency with file format
            >>    N of columns and input therein
            >>    any lines that do not contain noes
           
            >>    set appropriate class variables which will be used by NOE reading method
            '''            
            # here MAKE SURE THAT The input format holds >> last two columns must be Intensity[-2] and S/N[-1]
            #IF ENTIRE FORMAT CHECKED ABOVE, all should hold though!
            pass
        elif self.input_noe_format == 'vnmr':
            '''
            // TO DO
            self.ConvertToSparkyFormat('vnmr')
            '''
            pass
        elif self.input_noe_format == 'xeasy':
            '''
            // TO DO
            self.ConvertToSparkyFormat('xeasy')
            '''
            pass
        else:
            print 'ERROR! In input NMR files'
            print '>> Unrecognized or mistyped NOE list format provided, please choose one of the supported formats'
            sys.exit(0)      
        
    def SetGivenInput(self):
        
        #if self.input_peak_list:
            #self.ReadPeakList()
            
        #else:
            #print "No input peak list, simulated NOEs ... ?"        
        
        if self.input_noe:
            print "Reading from the NOE input file"
            noes = self.ReadInputNOEList()
            
        elif self.read_extracted_NOEs:
            print "Reading from a list of extracted NOEs"
            noes = self.ReadExtractedNOEs()
        else:
            print "Problem with NOE settings; check 'read_NOEs_from_file' and 'read_extracted_NOEs'; Only one of the options should be on!"
            sys.exit(0)
                        
        if self.check_recip and not self.check_intensity:
            print "Checking only for NOE reciprocity"
            noes = self.CheckRecipNOEs(noes)
        
        elif self.check_intensity and not self.check_recip:
            """
            check of intensity assumes the check of reciprocity as well
            """
            print "Checking for NOE intensity and (therefore also) reciprocity"
            self.check_recip = True
            noes = self.CheckRecipNOEs(noes)
            """
            //TO DO//
            Intensity check        
            """
        elif self.check_recip and self.check_intensity:
            print "Checking for NOE reciprocity and intensity"
            noes = self.CheckRecipNOEs(noes)
            """
            //TO DO//
            Intensity check        
            """
        else:        
            print "No NOE checking done..." 
    
        if self.merge:
            noes = self.MergeNOEs(noes)
        else:
            print 'Keeping ProR and ProS methyl groups of LEU/VAL as separate nodes ...'

        if self.merge_label_LV:
            print "No discrimination between the Leu / Val residue types ..."
            noes = self.MergeLabelLV(noes)
        else:
            print "Discriminating between Leu / Val residue types ..."
            
        if self.ligand_noes:
            self.lig_noes= self.ReadLigandNOEs()
        
        return noes     
    
    def SetNOEDict(self):
        self.noes = self.SetGivenInput()

    def ReadLigandNOEs(self):
        # read in 2 column file containing inter-molecular NOEs between protein and ligand
        # collect ligand noes to a dictionary
        try:
            fin = open(self.ligand_noe_file,"r")
        except:
            raise Exception('ERROR: file %s not found!'%self.ligand_noe_file)            
        lig_noes = {}
        for line in fin:
            stripped = line.strip()
            splitted = stripped.split()
            if len(splitted)==2:
                noe_id1 = splitted[0][1:]+splitted[0][0]
                noe_id2 = splitted[1][1:]+splitted[1][0]
                lig_noes.setdefault(noe_id1,[]).append(noe_id2)
                lig_noes.setdefault(noe_id2,[]).append(noe_id1)
            else:
                raise Exception("Format of the 'protein_ligand_noe_file' is incorrect. Two-column file expected, instead got % columns"%len(splitted))
        return lig_noes

    def MergeLabelLV(self,noe_dict):
        #print noe_dict
        merged_label_noes = {}
        for peak,conns in noe_dict.items():
            if peak[-1] == "V":
                relabelled_peak = peak[:-1]+"L"
            else:
                relabelled_peak = peak
            for c in conns:
                if c[-1] == "V":
                    relabelled_c = c[:-1]+"L"
                else:
                    relabelled_c = c
            merged_label_noes.setdefault(relabelled_peak,[]).append(relabelled_c)
        return merged_label_noes
    
    def MergeNOEs(self,noe_dict):
        # given the noes, merge those coming from the same residue but different methyl groups (ProR/ProS)
        print 'Merging NOEs from ProR/ProS methyl groups of the same Leu/Val residue'
        merged_noes = {}
        for peak,rnoes in noe_dict.items():
            try:
                new_peak = peak[:peak.index('_')]+peak[-1]
            except ValueError:
                new_peak = peak
            merge = []
            for noe in rnoes:
                try:
                    new_noe = noe[:noe.index('_')]+noe[-1]
                except ValueError:
                    new_noe = noe
                merge.append(new_noe)    
            merged_noes.setdefault(new_peak,merge)  #    reset merge
        return merged_noes
                        
    def CheckRecipNOEs(self,init_noes):
        #    check NOEs according to the flags given in the input file
        #    reciprocity means: NOE methyl1 >> methyl2; methyl2 >> methyl1 both exist
        #    only in case NOE is reciprocated it is kept;
        #    additional checks can be done    >>    intensity and signal to noise ratio (those assume reciprocity has been checked for)
        #    minimal check of the NOEs >> NOE is reciprocated 
        checked_noes = {}
        for peak,conns in init_noes.items():
            simple_conns = [c[0] for c in conns]    # check only for reciprocity, excluding the intensity and s/n check
            checked_noes.setdefault(peak,simple_conns)
        checked_noes = self.CheckReciprocity(checked_noes)
        return checked_noes
    
    def CheckReciprocity(self,dict_conns):
        # checks reciprocity of connection between 2 peaks
        # returns only connections that are reciprocated
                
        recip_conn_dict = {}
        all_values = list(set([v for value in dict_conns.values() for v in value]))
        indx_values = sorted(all_values)
        # create numpy array of string identifiers of peaks
        indx_values = np.array(indx_values)
        # create 2D array of zeros that contains all elements (peaks)
        conn_array = np.zeros(shape=(len(indx_values),len(indx_values)))
        # fill up 2D array of zeros to become an adjacency matrix where matrix[x,y] = 1 in case there is connection between x,y and 0 otherwise
        for key,conns in dict_conns.items():
            xindx = np.where(indx_values==key)
            for c in conns:
                yindx = np.where(indx_values==c)
                conn_array[xindx,yindx]=1   # conn array is filled up as an adjacency matrix
        
        recip = np.logical_and(conn_array,conn_array.T) # check where matrix and its transpose are both true (reciprocity) 
        recip_index = np.where(recip==True) # get the indices where both true
        
        recip1 = indx_values[recip_index[0]]    # get the values of NOEs1
        recip2 = indx_values[recip_index[1]]    # get the values of NOEs2 (reciprocals)
        
        for x in xrange(len(recip1)):
            recip_conn_dict.setdefault(recip1[x],[]).append(recip2[x])
        return recip_conn_dict


    def ReadExtractedNOEs(self):
        unfiltered_noes = {}
        fo = open(self.extracted_NOEs,'r')
        for line in fo:
            stripped = line.strip()
            splitted = stripped.split()
            if len(splitted) == 2:
                ID_int1 = re.findall(r'\d+',splitted[0])
                ID_type1 = splitted[0][0]
                newID_1 = str(ID_int1[0]) + ID_type1
                
                ID_int2 = re.findall(r'\d+',splitted[1])
                ID_type2 = splitted[1][0]
                newID_2 = str(ID_int2[0]) + ID_type2

                unfiltered_noes.setdefault(newID_1,[]).append((newID_2,None,None)) 
                unfiltered_noes.setdefault(newID_2,[]).append((newID_1,None,None))
            else:
                continue
            
        return unfiltered_noes
    
    def InputDimensions(self):
        #    set the order of dimensions according to the NOESY experiment
        if self.noesy_dimension == 3:            
            self.order_dimensions = 'C-C-H'
        elif self.noesy_dimension == 4:
            self.order_dimensions = 'C-H-C-H'
        else:
            raise Exception('Unrecognised dimensionality of the NOESY spectrum')
        
    def GetMethylType(self,peak):
        #    given peak in question, find its integer and str identifier
        #    returns None for str identifier if none there
        ID_int = re.findall(r'\d+', peak)
        for id in self.methyl_IDs:    
            try:
                indx = peak.index(id)
                ID_str = peak[indx]
                return ID_str,ID_int
            except Exception, e:
                continue
    
    def GetMethylID(self,numbers,type):
        if len(numbers) == 1:
            n= numbers[0]
            methyl = str(n) + type
        elif len(numbers) == 2:
            n1 = numbers[0]
            n2 = numbers[1]
            methyl = str(n1)+'_'+str(n2) + type
        elif len(numbers) == 0:
            raise Exception ('>> ERROR in peak input - no peak identifier encountered')
        elif len(numbers) > 2:
            raise Exception ('>> ERROR in peak input -- multiple integer peak identifiers encountered')
        else:
            raise Exception('>> ERROR in peak input - multiple peak identifiers encountered')

        return methyl
    
    def DetermineMethylID(self,freqC,freqH):
        #    take into account the methyl frequencies!
        '''// TO DO //'''
        pass

    def ReadPeakList(self):
        #    read in the peak list
        #    determine types of peaks
        methyl_list = open(self.input_peak_list,'r')
        for line in methyl_list:
                stripped = line.strip()
                splitted = stripped.split()
                if len(splitted) == 3:
                    peakID = splitted[0]
                    ID_str,ID_int = self.GetMethylType(peakID)
                    if ID_str == None:
                        print peakID
                        print '>> No methyl types given >> Implement the method DetermineMethylID!!'
                        print '>> Exiting ... '
                        
                        self.DetermineMethylID(splitted[1],splitted[2])
                    else:
                        methyl = self.GetMethylID(ID_int,ID_str)
                        self.peaks.append(methyl)
                else:
                    continue
        self.peaks = list(set(self.peaks))
        self.SortPeaks()

    def SortPeaks(self):
        # sort the list of peaks according to the first integer identifier
        reduced_list = [int(re.findall(r'\d+', peak)[0]) for peak in self.peaks]
        self.peaks = [x for (y,x) in sorted(zip(reduced_list,self.peaks))]    
    
    def ReadInputNOEList(self):
        # first sets the expected order of the dimensions (given in the input file)
        # reads in an NOE list 
        # determine types of peaks
        # checks for dimensionality and consistency with file format should have been done by this point!
        # expects the noe identifiers to be given in the first column and to match the input dimension 

        unfiltered_noes = {}
        print 'Reading NOE list ... '
        self.InputDimensions()  
        frequencies_order = self.order_dimensions.split('-')   
        noe = open(self.input_noe_list,'r')
        for line in noe:
            stripped = line.strip()
            splitted = stripped.split()
            noe_label = splitted[0].split('-')
            # if auto peak, skip
            if noe_label[-1] == frequencies_order[-1] and noe_label[-2] == frequencies_order[-2]:
                continue
            
            if len(frequencies_order) == 3:
                noe1_str,noe1_int= self.GetMethylType(noe_label[0])
                noe2_str,noe2_int= self.GetMethylType(noe_label[1])
            
                '''
                if noe1_str == None:
                    self.DetermineMethylID(splitted[1],splitted[3])
                if noe2_str == None:
                    self.DetermineMethylID(splitted[2])
                    
                >> depends on the type of the experiment --> would be good to do this only for
                the HMQC peaks and then relate those to the NOEs
                '''
            elif len(frequencies_order) == 4:
                noe1_str,noe1_int= self.GetMethylType(noe_label[0])
                noe2_str,noe2_int= self.GetMethylType(noe_label[2])
                   #print noe_label                
                '''
                if noe1_str == None:
                    self.DetermineMethylID(splitted[1],splitted[2])
                if noe2_str == None:
                    self.DetermineMethylID(splitted[3],splitted[4])
                    
                >> depends on the type of the experiment --> probably would be good to do this only for
                the HMQC peaks and then relate those to the NOEs
                '''
            else:
                print 'Incorrect or unrecognised dimensionality of the NOESY spectrum provided!' 
            
            methyl1 = self.GetMethylID(noe1_int,noe1_str)
            methyl2 = self.GetMethylID(noe2_int,noe2_str)
            
            #print methyl1, ' >> ',methyl2                                           
            #unfiltered_noes.setdefault(methyl1,[]).append((methyl2,float(splitted[-2]),float(splitted[-1])))
            unfiltered_noes.setdefault(methyl1,[]).append((methyl2,float(0),float(0)))

        return unfiltered_noes

    def CheckInt_StoN(self):
        '''
        // TO DO//
        '''
        pass
            
    def AssembleResults(self,resultdict,short,long,assign):
        # result dict is nested dictionary
        try:
            resultdict[long]
        except KeyError:
            resultdict.setdefault(long,{})
        try:
            resultdict[long][short]
            for peak,atom in assign.items():
                resultdict[long][short].setdefault(peak,[]).append(atom)
        except KeyError:
            resultdict[long].setdefault(short,{})
            for peak,atom in assign.items():
                resultdict[long][short].setdefault(peak,[]).append(atom)                    
        return resultdict       
                
    def SetEdgeAssignment(self,medges_update,assign,snodes,sedges,bnodes,bedges):
        for snode,bnode in assign.items():
            e1,e2= self.GetIntEdges(snode,bnode,snodes,bnodes, sedges,bedges)
            medges_update= self.EdgeAssignment(medges_update,e1,e2)
        escore = self.EdgeScore(medges_update)
        return escore
        
    def EdgeAssignment(self,medges,edges1,edges2):        
        for edge1 in edges1:
            for col in range(len(medges[edge1])):
                if medges[edge1,col]==1 and col not in edges2:   # if there is an edge for the node in G1 and there is no edge for the corresponding node in G2
                    medges[edge1,col] = 0  # set this correspondence in edges to 0 (i.e. in the current mapping edge1 of G1 can not correspond to some of the edges in 
        return medges
    
    def EdgeScore(self,medges):
        nonzero = 0 
        for row in medges:
            if np.count_nonzero(row) != 0:  # if entire row is not set to 0
                nonzero+=1  # the number of self.edgesleft is equal to the number of rows with non zero elements
        return nonzero
                

    
