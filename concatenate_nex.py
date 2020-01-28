from Bio.Nexus import Nexus
import sys

#usage : python concatenate_nex.py filelist.txt

filename = sys.argv[1] # provide a list of files in nexus format as text file as command line argument
with open(filename, 'r') as f:
    file_list = f.read().split('\n')
for i in file_list: # print the list of input files 
    print(i)

nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]
combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open('concatenated.nex', 'w')) # writes output file named as concatenated.nex