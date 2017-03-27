MAGMA provides a collection of methods to define simple, undirected, labelled graphs on the basis of input structure (PDB) files and NMR data.

MAGMA main features:
* importing of PDB files, possibly containing multiple conformations (e.g. multi PDB)
* extraction of relevant atomic coordinates and generation of pseudoatoms (when so specified)
* generation of graph objects on the basis of atomic coordinates and inter-atomic distances
* generation of graph objects on the basis of chemical shift values and nuclear Overhauser enhancement rates' chemical shift values
* traversal of graph vertices using depth- or breadth-first search
* evaluation of similarity between the graphs using subgraph_isomorphism algorithm and/or MCES algorithm
* generation of mapping from the vertices of one graph (i.e. nmr data graph) to the vertices of another graph (i.e. structure graph) on the basis of subgraph isomorphism and/or maximal common edge subgraph

MAGMA program, documentation and examples [here](http://magma.chem.ox.ac.uk).


## RUNNING MAGMA

MAGMA can be run using Linux/Unix and Windows.

For a successful operation of MAGMA, you need to have:
		
* GNU Compiler Collection (GCC)
* python (__version__= '2.7' or higher)
		
 And the following python packages installed on your system.
 
* numpy (__version__= '1.8.1' or higher)
* matplotlib (__version__ = '1.3.1' or higher)
* igraph (__version__= '0.7.1' or higher)
* networkx (__version__= '1.8.1' or higher)
* munkres (__version__= '1.0.7' or higher)
			
MAGMA python-version is run from a command line as follows:
 `python execute_magma.py /path/to/input/text/file'

## MAGMA OUTPUTS

MAGMA calculation will create a subdirectory called `results` in the directory where `input_protein_name.txt` file is located

After a successful calculation the directory `results` will contain:
								
* `combined_result.res`

This is the FINAL assignment file.
In the first column 1H-13C peak identifier is given.
In the second and any nth column following ":" all residue assignment options for the peak are given

The rest of the outputs contain details of certain steps of the overall MAGMA protocol.
These are less relevant for the user, but useful for the analysis of run progression.

* `hungarian_assignmentX.out`, the assignment guess prior to execution of the exact MCES algorithm (relevant to the user!). In case the calculation takes very long to converge -- this assignment can be treated as the best guess, approximate assignment.  X - full for "all_subgraphs" mode. X - subgraph numerical id for "connected_subgraphs" mode
															
* `conv.2.x.out`, text file monitoring details of the first step of MAGMA optimization prior to execution of complete MCES search. For each step of iterative optimization (column 1) enumerate how many vertices (column 3) have scored a particular score (column 2)

* `conv.1.x.out`. text file monitoring progress of the second step of MAGMA optimization prior to execution of complete MCES search. For each step of iterative optimization (column 1) store the score of collected mces (column 2)
			
* `mces_n.txt`, lists all mces solutions found by McGregor algorithm for the nth connected subgraph of the nmr data graph

* `mces_n.txt.B`

* `mces_n.txt.G`

* `vf2_n.txt`, consolidation of assignment through subgraph isomorphism, using VF2 algorithm, for the nth connected subgraph of the nmr data graph


## CREDITS

author: Iva Pritisanac (iva.pritisanac[at]gmail.com), University of Oxford

other contributions:

* classes Parser, Molecule and Structure -- Matteo Degiacomi (matteo.degiacomi[at]gmail.com), University of Oxford
* functions for monitoring optimization progress and producing results summary files -- Andrew Baldwin (andrew.baldwin[at]chem.ox.ac.uk), University of Oxford
* implementation of subgraph isomorphism algorithm from [igraph](igraph.org)
* implementation of Hungarian algorithm from the python package [munkres 1.0.7](https://pypi.python.org/pypi/munkres/)
* graph objects definitions, manipulations, and plotting python package [networkx 1.9.1](https://pypi.python.org/pypi/networkx/)
* implementation of algorithm for graph breadth first search python package [networkx 1.9.1](https://pypi.python.org/pypi/networkx/)
  

## TO DO/ONGOING DEVELOPMENT

* automatic iteration over a user specified range of distance thresholds to define an optimal threshold. Status: ONGOING.
* a call to mcgregor method to enable mode merge_proRS 'off' mode after the calculation in mode merge_proRS 'on'. Status: ONGOING
* cython the slowest python functions 

## KNOWN BUGS

* method remove_small_subgraphs of class Graph occasionally fails (unknown reason). Current treatment: Exception
* method reprioritize of class McGregor (unknown reason). Current treatment: None


## NOTE TO USERS

To report bugs or suggest developments please contact: iva.pritisanac[at]gmail.com
Have fun using MAGMA :-)