import os,sys


NOE_input_path="/home/beast/ipritisanac/Benchmark_input_files/msrB/filtered_merged_NOEs.txt"
base_input_path="/home/beast/ipritisanac/Benchmark_input_files/msrB/"
N_prior=10
N_mces=3

dir = "joined_calc_input/"
mces_mode = "one"
mode = "all_subgraphs"

for i in os.listdir(dir):
    fin = open(dir+i,"r")
    fout=open(dir+i[:-5]+"jo.txt","w")
    for line in fin:
        splitted = line.split()
	if splitted[0]=="extracted_NOE_list":
	    fout.write("%s\t%s\n"%("extracted_NOE_list",NOE_input_path))
        elif splitted[0]=="output_file_McGregor":
            fout.write("%s\t%s\n"%("output_file_McGregor",base_input_path+"joined_results/mces_"+i[:-5]))
	elif splitted[0]=="output_file_vf2":
            fout.write("%s\t%s\n"%("output_file_vf2",base_input_path+"joined_results/vf2_"+i[:-5]))
        elif splitted[0]=="priority_iterations":
            fout.write("%s\t%s\n"%("priority_iterations",str(N_prior)))
        elif splitted[0]=="N_priority_MCES":
            fout.write("%s\t%s\n"%("N_priority_MCES",str(N_mces)))
	elif splitted[0]=="mces_mode":
            fout.write("%s\t%s\n"%("mces_mode",mces_mode))
	elif splitted[0]=="mode":
            fout.write("%s\t%s\n"%("mode",str(mode)))
        else:
            fout.write(line)
    fin.close()
    fout.close()
