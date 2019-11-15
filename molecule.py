"""
__author__: Matteo Degiacomi
PDB file loading and manipulation. Extracted from BioBOx, (c) Matteo Degiacomi
"""

import numpy as np
import scipy.spatial.distance as S
from copy import deepcopy

## Read and manipulate a PDB file.
class Molecule():

    ## initialize Molecule properties data structure.
    # Store properties associated to protein in properties dictionary entry "data". Every information about every atom
    # is stored in this order (as string): [atom/hetatm, index, atomname, resname, chain name, residue ID, beta factor, occupancy, atomtype].
    # self.knowledge contains a knowledge base about atoms and residues properties. Default values are:
    # - 'residue_mass' property stores the average mass for most common aminoacids (values from Expasy website)
    # - 'atom_vdw' vdw radius of common atoms
    # - 'atom_mass' mass of common atoms
    # the knowledge base can be edited. For instance, to add information about residue "TST" mass in molecule M type: M.knowledge['mass_residue']["TST"]=142.42
    def __init__(self,p=np.array([[],[]]),r=1.9):

        #super(Molecule,self).__init__(r=np.array([]))
        if p.ndim==3:
            ##self.coordinates numpy array containing an ensemble of alternative coordinates in 3D space
            self.coordinates=p
        elif p.ndim==2:
            self.coordinates=np.array([p])
        else:
            print("ERROR: expected numpy array with 2 or three dimensions, but %s dimensions were found"%p.ndim)
            return -1

        ##index of currently selected conformation.
        self.current=0
        ##pointer to currently selected conformation.
        self.points=self.coordinates[self.current]


        ##collection of properties. By default, 'center' (geometric center of the Structure) is defined
        self.properties={}
        #self.properties['center']=self.get_center()
        self.properties['radius']=r #average radius of every point



        #pdb data information. for every atom store (in order): ATOM/HETATM, index, atomname, resname, chain name, residue ID, occupancy, beta factor, atomtype
        self.properties['data']=np.array([])

        ##knowledge base about atoms and residues properties (entry keys: 'residue_mass', 'atom_vdw', 'atom_,mass' can be edited)
        self.knowledge={}
        self.knowledge['atom_vdw']={'H': 1.20, 'N': 1.55, 'NA': 2.27, 'CU': 1.40, 'CL': 1.75, 'C': 1.70, 'O': 1.52, 'I': 1.98, 'P': 1.80, 'B': 1.85, 'BR': 1.85, 'S': 1.80, 'SE': 1.90,\
                                 'F': 1.47, 'FE': 1.80, 'K':  2.75, 'MN': 1.73, 'MG': 1.73, 'ZN': 1.39, 'HG': 1.8, 'XE': 1.8, 'AU': 1.8, 'LI': 1.8, '.': 1.8}

    
    
    ##return property from properties array
    #@param prop desired property to extract from property dictionary
    #@retval or nan if failed
    def get(self,prop):
        if str(prop) in self.properties:
            return self.properties[str(prop)]
        else:
            print("property %s not found!"%prop)
            return float('nan')

    ##select current frame (place frame pointer at desired position)
    # @param pos number of alternative conformation (starting from 0)
    def set_current(self, pos):
        if pos<self.coordinates.shape[0]:
            self.current=pos
            self.points=self.coordinates[self.current]
            #self.properties['center']=self.get_center()
        else:
            print("ERROR: position %s requested, but only %s conformations available"%(pos,self.coordinates.shape[0]))

    ##create a new property.
    def add_property(self,name,value):
        self.properties[name]=value
    
    ## get points coordinates.
    # @param indices indices of points to select. If none is provided, all points coordinates are returned.
    # @retval coordinates of all points indexed by the provided indices list, or all of them if no list is provided.
    def get_xyz(self, indices=[]):
        if indices==[]:
                return self.points
        else:
                return self.points[indices]

    ## set point coordinates.
    # @param coords array of 3D points
    def set_xyz(self,coords):
        self.coordinates[self.current]=deepcopy(coords)
        self.points=self.coordinates[self.current]

    ## add a new alternative conformation to the database
    # @param coords array of 3D points, or array of arrays of 3D points (in case multiple alternative coordinates must be added at the same time)
    def add_xyz(self, coords):
        #self.coordinates numpy array containing an ensemble of alternative coordinates in 3D space

        if self.coordinates.size==0 and coords.ndim==3:
                self.coordinates=deepcopy(coords)
                self.set_current(0)

        elif self.coordinates.size==0 and coords.ndim==2:
                self.coordinates=deepcopy(np.array([coords]))
                self.set_current(0)

        elif self.coordinates.size>0 and coords.ndim==3:
                self.coordinates=np.concatenate((self.coordinates, coords))
                self.set_current(self.current+1) # set new frame to the first of the newly inserted ones

        elif self.coordinates.size>0 and coords.ndim==2:
                self.coordinates=np.concatenate((self.coordinates, np.array([coords])))
                self.set_current(self.current+1) # set new frame to the first of the newly inserted one

        else:
                print("ERROR: expected numpy array with 2 or three dimensions, but %s dimensions were found"%np.ndim)
                return -1    
    
    ##return information from knowledge base
    #@param prop desired property to extract from knowledge base
    #@retval value associated to requested propery, or nan if failed
    def know(self,prop):
        if str(prop) in self.knowledge:
                return self.knowledge[str(prop)]
        else:
                print("entry %s not found in knowledge base!"%prop)
                return float('nan')


    ## read a pdb (possibly containing containing multiple models).
    # models are split according to ENDMDL and END statement. All alternative coordinates are expected to have the same atoms. After loading, the first model (M.current_model=0) will be set as active.
    # @param pdb PDB filename
    # @param include_hetatm if True, HETATM will be included (they get skipped if False)
    def import_pdb(self,pdb,include_hetatm=False):

        try:
            f_in=open(pdb,"r")
        except Exception as E:
            e=E
            print("%s"%e)
            raise Exception('ERROR: file %s not found!'%pdb)

        #store filename
        self.add_property("filename",pdb)

        data_in=[]
        p=[]
        r=[]
        c=[]
        alternative=[]
        biomt=[]
        for line in f_in:
            record=line[0:6].strip()

            
            #load symmetry matrix, if any is present
            if "REMARK 350   BIOMT" in line:
                try:
                    biomt.append(line.split()[4:8])
                except:
                    raise Exception("ERROR: biomatrix format seems corrupted")
            

            #if a complete model was parsed store all the saved data into self.properties entries (if needed) and temporary alternative coordinates list
            if record=="ENDMDL" or record=="END":

                if len(alternative)==0:

                        #load all the parsed data in superclass properties['data'] and points data structures
                        try:
                            #ATOM/HETATM, index, atomname, resname, chain name, residue ID, beta factor, occupancy, atomtype
                            self.properties['data']=np.array(data_in).astype(str)
                        except:
                            raise Exception('ERROR: something went wrong when loading the structure %s!\nERROR: are all the columns separated?'%pdb)

                        #saving vdw radii	
                        try:
                            self.properties['radius']=np.array(r)
                        except:
                            raise Exception('ERROR: something went wrong when loading the structure %s!\nERROR: are all the columns separated?'%pdb)

                #save 3D coordinates of every atom and restart the accumulator
                try:
                    if len(p)>0:
                        alternative.append(np.array(p))
                    p=[]
                except:
                    raise Exception('ERROR: something went wrong when loading the structure %s!\nERROR: are all the columns separated?'%pdb)


            if record=='ATOM' or (include_hetatm and record=='HETATM'):

                #extract xyz coordinates (save in list of point coordinates)
                p.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])

                #if no complete model has been yet parsed, load also information about atoms(resid, resname,...)
                if len(alternative)==0:
                        w=[]
                        w.append(line[0:6].strip()) #extract ATOM/HETATM statement
                        w.append(line[6:11].strip()) #extract atom index
                        w.append(line[12:17].strip()) #extract atomname
                        w.append(line[17:20].strip()) #extract resname
                        w.append(line[21].strip()) #extract chain name
                        w.append(line[22:26].strip()) #extract residue id

                        #extract occupancy
                        try:
                            w.append(float(line[54:60]))
                        except:
                            w.append(1.0)

                        #extract beta factor
                        try:
                            #w.append("{0.2f}".format(float(line[60:66])))
                            w.append(float(line[60:66]))
                        except:
                            w.append(0.0)

                        #extract atomtype
                        try:
                            w.append(line[76:78].strip())
                        except:
                            w.append("")

                        #use atomtype to extract vdw radius
                        try:
                            r.append(self.know('atom_vdw')[line[76:78].strip()])
                        except:
                            r.append(self.know('atom_vdw')['.'])

                        data_in.append(w)

        f_in.close()




        #if p list is not empty, that means that the PDB file does not finish with an END statement (like the ones generated by SBT, for instance).
        #In this case, dump all the remaining stuff into alternate coordinates array and (if needed) into properties dictionary.
        if len(p)>0:

                #if no model has been yet loaded, save also information in properties dictionary.
                if len(alternative)==0:

                        #load all the parsed data in superclass properties['data'] and points data structures
                        try:
                            #ATOM/HETATM, index, atomname, resname, chain name, residue ID, beta factor, occupancy, atomtype
                            self.properties['data']=np.array(data_in).astype(str)
                        except:
                            raise Exception('ERROR: something went wrong when saving data in %s!\nERROR: are all the columns separated?'%pdb)

                        #saving vdw radii
                        try:
                            self.properties['radius']=np.array(r)
                        except:
                            raise Exception('ERROR: something went wrong when saving van der Waals radii in %s!\nERROR: are all the columns separated?'%pdb)


                #save 3D coordinates of every atom and restart the accumulator
                try:
                    if len(p)>0:
                            alternative.append(np.array(p))
                    p=[]
                except:
                    raise Exception('ERROR: something went wrong when saving coordinates in %s!\nERROR: are all the columns separated?'%pdb)
        
        #transform the alternative temporary list into a nice multiple coordinates array
        if len(alternative)>0:
                try:
                    alternative_xyz=np.array(alternative).astype(float)
                except:
                    alternative_xyz=np.array([alternative[0]]).astype(float)
                    print('WARNING: found %s models, but their atom count differs'%len(alternative))
                    print('WARNING: treating only the first model in file %s'%pdb)
                    #raise Exception('ERROR: models appear not to have the same amount of atoms')

                self.add_xyz(alternative_xyz)

        else:
                raise Exception('ERROR: something went wrong when saving alternative coordinates in %s!\nERROR: no model was loaded... are ENDMDL statements there?'%pdb)


    ## select specific atoms in the protein providing chain, residue ID and atom name.
    # @param chain selection of a specific chain name (accepts * as wildcard). Can also be a list or numpy array of strings.
    # @param res residue ID of desired atoms (accepts * as wildcard). Can also be a list or numpy array of of int.
    # @param atom name of desired atom (accepts * as wildcard). Can also be a list or numpy array of strings.
    # @param get_index if set to True, returns the indices of selected atoms in self.points array (and self.properties['data'])
    # @param use_resname if set to True, consider information in "res" variable as resnames, and not resids
    # @retval coordinates of the selected points and, if get_index is set to true, their indices in self.points array.
    def atomselect(self,chain,res,atom,get_index=False,use_resname=False):

        #chain name boolean selector
        if type(chain) is str:
            if chain=='*':
                chain_query=np.array([True]*len(self.points))
            else:
                chain_query=self.properties['data'][:,4]==chain

        elif type(chain) is list or type(chain).__module__=='numpy':
            chain_query=self.properties['data'][:,4]==chain[0]
            for c in range(1,len(chain),1):
                chain_query=np.logical_or(chain_query,self.properties['data'][:,4]==chain[c])
        else:
            print("ERROR: wrong type for chain selection. Should be str, list, or numpy")
            return -1


        #resid boolean selector
        #pick appropriate column for boolean selection
        if use_resname:
            refcolumn=3
        else:
            refcolumn=5

        if type(res) is str:
            if res=='*':
                res_query=np.array([True]*len(self.points))
            else:
                res_query=self.properties['data'][:,refcolumn]==res
        elif type(res) is int:
            res_query=self.properties['data'][:,refcolumn]==str(res)
        elif type(res) is list or type(res).__module__=='numpy':
            res_query=self.properties['data'][:,refcolumn]==str(res[0])
            for r in range(1,len(res),1):
                res_query=np.logical_or(res_query,self.properties['data'][:,refcolumn]==str(res[r]))   
        else:
            print("ERROR: wrong type for resid selection. Should be int, list, or numpy")
            return -1

        #atom name boolean selector
        if type(atom) is str:
            if atom=='*':
                atom_query=np.array([True]*len(self.points))
            else:  
                atom_query=self.properties['data'][:,2]==atom
        elif type(atom) is list or type(atom).__module__=='numpy':
            atom_query=self.properties['data'][:,2]==atom[0]
            for a in range(1,len(atom),1):
                atom_query=np.logical_or(atom_query,self.properties['data'][:,2]==atom[a])   
        else:
            print("ERROR: wrong type for atom selection. Should be str, list, or numpy")
            return -1    

        #slice data array and return result (colums 5 to 7 contain xyz coords)
        query=np.logical_and(np.logical_and(chain_query,res_query),atom_query)

        #@todo mass of termini to add!
        if get_index==True:
            return [self.points[query],np.where(query==True)[0]]
        else:
            return self.points[query]



    ##select atoms having the same residue and chain as a given atom (or list of atoms)
    # @param index indices of atoms of choice (integer of list of integers)
    # @param get_index if set to True, returns the indices of selected atoms in self.points array (and self.properties['data'])
    # @retval coordinates of the selected points and, if get_index is set to true, their indices in self.points array.
    def same_residue(self,index,get_index=False):

        D=self.properties['data']
        l=D[index]

        if len(l.shape)==1:
            l=l.reshape(1,len(l))

        test=np.logical_and(D[:,4]==l[:,4],D[:,5]==l[:,5]) 

        idxs=np.where(test)[0]
        if len(idxs)>0:
            pts=self.points[idxs]
        else:
            pts=[]

        if get_index==True:
            return pts,idxs
        else:
            return pts



    ## aggregate data and point coordinates, and return in a unique data structure (list, str for points data, float for point coordinates)
    # @retval list aggregated data and coordinates for every point, as string. Same order as a pdb file, i.e. ATOM/HETATM, index, atomname, resname, chain name, residue ID, x, y, z, occupancy, beta factor, atomtype
    def get_pdb_data(self,index=""):
        
        if index=="":
                index=range(0,len(self.points),1)

        #create a list containing all infos contained in pdb (point coordinates and properties)
        d=[]
        for i in index:
            d.append([self.properties['data'][i,0],self.properties['data'][i,1],self.properties['data'][i,2],self.properties['data'][i,3],self.properties['data'][i,4],self.properties['data'][i,5],self.points[i,0],self.points[i,1],self.points[i,2],self.properties['data'][i,6],self.properties['data'][i,7],self.properties['data'][i,8]])
        
        return d


    ## overload superclass method for writing (multi)pdb.
    # @param outname name of pdb file to be generated.
    # @param index indices of atoms to write to file. If empty, all atoms are returned. Index values obtaineable with a call like: index=molecule.atomselect("A",[1,2,3],"CA",True)[1]
    # @param conformations list of conformation indices to write to file. By default, a multipdb with all conformations will be produced.
    def write_pdb(self,outname, conformations=[], index=""):

        #store current frame, so it will be reestablished after file output is complete
        currentbkp=self.current

        #if a subset of all available frames is requested to be written, select them first
        if len(conformations)==0:
                frames=range(0,len(self.coordinates),1)
        else:
                if np.max(conformations)<len(self.coordinates):
                        frames=conformations
                else:
                        print("ERROR: requested coordinate index %s, but only %s are available"%(np.max(conformations), len(self.coordinates)))
                        return -1

        f_out=open(outname,"w")

        for f in frames:

                #get all informations from PDB (for current conformation) in a list
                self.set_current(f)
                d=self.get_pdb_data(index)

                for i in range(0,len(d),1):
                    #create and write PDB line
                    L='%-6s%5s  %-4s%-4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%(d[i][0],d[i][1],d[i][2],d[i][3],d[i][4],d[i][5],float(d[i][6]),float(d[i][7]),float(d[i][8]),float(d[i][9]),float(d[i][10]),d[i][11])
                    f_out.write(L)

                f_out.write("END\n")

        f_out.close()

        self.set_current(currentbkp)

        return


    ##given a list of indices, compute the all-vs-all distance and return only couples below a given cutoff distance
    #useful for the detection of disulfide bridges or linkable sites via cross-linking (approximation, supposing euclidean distances)'
    #@param idx indices of atoms to check.
    #@param cutoff minimal distance to consider a couple as linkable.
    #@retval nx3 numpy array containing, for every valid connection, id of first atom, id of second atom and distance between the two.
    def get_couples(self, idx, cutoff):

        points1=self.get_xyz()[idx]
        dist=S.cdist(points1,points1)
        couples=np.array(np.where(dist<cutoff))

        res=[]
        for c in couples.transpose():
            if c[0]>c[1]:
                res.append([idx[c[0]], idx[c[1]], dist[c[0],c[1]]])

        return np.array(res)


