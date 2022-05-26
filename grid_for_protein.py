#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 22:50:02 2022

@author: dozeduck
"""
import getopt
import sys
# import string
# import numpy


args=sys.argv[1:]                                                                   # get arguments from usr's input; filename is sys.arg[0], so args start from [1:]


try:
   opts, args = getopt.getopt(args,"h:p:l:o:d:",["help",
                                                   "input_pdb=",
                                                   "input_lig=",
                                                   "number_per_sublibrary=",
                                                   "output_name=",
                                                   "distance_box="])
except getopt.GetoptError:
   print ('grid_for_protein.py -p <inputprotein> -l <inputligands> -o <output_name> -d <distance_between_box_protein>')
   sys.exit(2)                                                                      # Exiting the program raises a SystemExit exception, 0 means normal exit, and others are abnormal exits.
 
for opt, arg in opts:                                                               # Generate several pairs of value, e.g: opr,arg = -i,PLDXPAL
   if opt == '-h':
      print ('grid_for_protein.py -p <inputprotein> -l <inputligands> -o <output_name> -d <distance_between_box_protein>')
      sys.exit()
   elif opt in ("-p", "--input_pdb"):
      f = str(arg)
   elif opt in ("-l", "--input_lig"):
      l = str(arg)
   elif opt in ("-o", "--outputfile"):
      out_put_name = str(arg)
   elif opt in ("-d", "--distance"):
      distance = int(arg)

# Read Protein PDB and Protein size
X_peratom = []
Y_peratom = []
Z_peratom = []
file = open(f,'r')
for line in file:
    if(line.split()[0] == 'ATOM' or line.split()[0] == 'HETATOM'):
        X_peratom.append(float(line.split()[6]))
        Y_peratom.append(float(line.split()[7]))
        Z_peratom.append(float(line.split()[8]))

min_X = min(X_peratom)
max_X = max(X_peratom)
min_Y = min(Y_peratom)
max_Y = max(Y_peratom)
min_Z = min(Z_peratom)
max_Z = max(Z_peratom)
file.close()

# Read Ligands SDF and find the long axis
x_peratom = []
y_peratom = []
z_peratom = []
s = list(map(chr, range(ord('a'), ord('z') + 1)))                                   # generate alphabet
s1 = [] 
num_mol = 0                                                                            # capitalize all alphabet
for x in s:
    s1.append(x.upper())
file = open(l,'r')
for line in file:
    try:
        if(int(line.split()[5-15]) == 0):                                                
            x_peratom.append(float(line.split()[0]))
            y_peratom.append(float(line.split()[1]))
            z_peratom.append(float(line.split()[2]))
    except:
        pass
file.close()
file = open(l,'r')
for line in file: 
    try:                                                                 # count the number of molecules in this library
        if(line.split()[0] == "$$$$"):
            num_mol += 1
    except:
        pass
file.close()
# print(num_mol)
min_x = min(x_peratom)
max_x = max(x_peratom)
min_y = min(y_peratom)
max_y = max(y_peratom)
min_z = min(z_peratom)
max_z = max(z_peratom)

box_8_point = []                                                                    # generate the box for protein 
for x in [min_X, max_X]:
    for y in [min_Y, max_Y]:
        for z in [min_Z, max_Z]:
           box_8_point.append([x,y,z])                                              # add 8 points coord of the box
with open('box.pdb',"w") as box:
    box.write(str("HEADER     CORNERS  OF  BOX \n"))
    box.write(str("REMARK     CENTER (X Y Z)" + str(round(sum([min_X,max_X])/2,3)) + " " + str(round(sum([min_Y,max_Y])/2,3)) + " " + str(round(sum([min_Z,max_Z])/2,3)) + "\n"))
    box.write(str("REMARK     DIMENSIONS (X Y Z)" + " " + str(round(abs(min_X + max_X),3)) + " " + str(round(abs(min_Y + max_Y),3)) + " " + str(round(abs(min_Z + max_Z),3)) + "\n"))
    print("%4s%7d%5s%4s%6d%12.3f%8.3f%8.3f" %  ("ATOM" ,     # Formatted output, %4s, right-aligned, the output occupies 4 columns in total. If the length is less than 4 columns, the left end will be filled with spaces. If it is greater than 4 columns, the actual length will be output as a string
                                         1,                          # %7d, right-aligned, the output occupies a total of 7 columns, if the length is less than 7 columns, the left end is filled with spaces, signed decimal certificate integer
                                         "DU"+s1[0],                           # %-4s, left-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the right end is filled with spaces, if it is greater than 4 columns, the actual length is output as a string
                                         "BOX",                          # %1s, right-aligned, the output occupies a total of 1 column. If it is less than 1 column, it will be filled with spaces from the left end. If it is greater than 1 column, the actual length will be output as a string
                                         1,                            # %2s, right-aligned, the output occupies 2 columns in total. If it is less than 2 columns, it will be filled with spaces from the left end. If it is greater than 2 columns, the actual length will be output as a string
                                         box_8_point[0][0]-distance,                         # %4d, right-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the left end is filled with spaces, a signed decimal certificate integer
                                         box_8_point[0][1]-distance,                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                         box_8_point[0][2]-distance), file = box )
    print("%4s%7d%5s%4s%6d%12.3f%8.3f%8.3f" %  ("ATOM" , 2, "DU"+s1[1], "BOX",  1, box_8_point[1][0]-distance, box_8_point[1][1]-distance, box_8_point[1][2]+distance), file = box )
    print("%4s%7d%5s%4s%6d%12.3f%8.3f%8.3f" %  ("ATOM" , 3, "DU"+s1[2], "BOX",  1, box_8_point[2][0]-distance, box_8_point[2][1]+distance, box_8_point[2][2]-distance), file = box )
    print("%4s%7d%5s%4s%6d%12.3f%8.3f%8.3f" %  ("ATOM" , 4, "DU"+s1[3], "BOX",  1, box_8_point[3][0]-distance, box_8_point[3][1]+distance, box_8_point[3][2]+distance), file = box )
    print("%4s%7d%5s%4s%6d%12.3f%8.3f%8.3f" %  ("ATOM" , 5, "DU"+s1[4], "BOX",  1, box_8_point[4][0]+distance, box_8_point[4][1]-distance, box_8_point[4][2]-distance), file = box )
    print("%4s%7d%5s%4s%6d%12.3f%8.3f%8.3f" %  ("ATOM" , 6, "DU"+s1[5], "BOX",  1, box_8_point[5][0]+distance, box_8_point[5][1]-distance, box_8_point[5][2]+distance), file = box )
    print("%4s%7d%5s%4s%6d%12.3f%8.3f%8.3f" %  ("ATOM" , 7, "DU"+s1[6], "BOX",  1, box_8_point[6][0]+distance, box_8_point[6][1]+distance, box_8_point[6][2]-distance), file = box )
    print("%4s%7d%5s%4s%6d%12.3f%8.3f%8.3f" %  ("ATOM" , 8, "DU"+s1[7], "BOX",  1, box_8_point[7][0]+distance, box_8_point[7][1]+distance, box_8_point[7][2]+distance), file = box )

    print("%-6s%5d%5d%5d%5d" % ("CONECT",1,2,3,5), file = box)
    print("%-6s%5d%5d%5d%5d" % ("CONECT",2,1,4,6), file = box)
    print("%-6s%5d%5d%5d%5d" % ("CONECT",3,1,4,7), file = box)
    print("%-6s%5d%5d%5d%5d" % ("CONECT",4,2,3,8), file = box)
    print("%-6s%5d%5d%5d%5d" % ("CONECT",5,1,6,7), file = box)
    print("%-6s%5d%5d%5d%5d" % ("CONECT",6,2,5,8), file = box)
    print("%-6s%5d%5d%5d%5d" % ("CONECT",7,3,5,8), file = box)
    print("%-6s%5d%5d%5d%5d" % ("CONECT",8,4,6,7), file = box)
        


# write the average coordination for each molecules
file = open(l,'r')
count = 0
atom_amount_line = []
atom_amount = []                                                                     # include the number of atoms for each molecules
ligand_library = []
for line in file:
    ligand_library.append(line)
    try:
        if(line.split()[0] == "SeeSAR"):                                                
            atom_amount_line.append(count+1)
    except:
        pass
    count += 1
file.close()
for i in atom_amount_line:
   atom_amount.append(ligand_library[i].split()[0])
# average_coord = []

file = open(l,'r')                                                                                # adding binding affinity information
count = 0
low_affinity_amount_line = []
low_affinity_amount = [] 
high_affinity_amount_line = []                                                                    # include the number of atoms for each molecules
high_affinity_amount = []
which_compound = []  
number_of_molecules = 0                                                                             # save which compounds have the affinity
for line in file:
    try:
        if(line.split()[0] == "SeeSAR"):
            number_of_molecules += 1
    except:
        pass
    try:
        if(line.split()[1] == "<BIOSOLVEIT.HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY"): 
                                                          
            low_affinity_amount_line.append(count+1)
            high_affinity_amount_line.append(count+4)
            # print(ligand_library[low_affinity_amount_line[0]-1].split()[:]) 
            which_compound.append(number_of_molecules-1)   
    except:
        pass
    count += 1
which_compound.append(0)
file.close()

for i in low_affinity_amount_line:
    low_affinity_amount.append(float(ligand_library[i].split()[0]))
for i in high_affinity_amount_line:
    high_affinity_amount.append(float(ligand_library[i].split()[0]))
# print(low_affinity_amount[2]," ", high_affinity_amount[2])
# print(which_compound)
count = 0
b = 0
with open(out_put_name,"w") as dummy:
    for i in range(num_mol):
        if(i == which_compound[b]):
        # print(count,count+int(atom_amount[i]))
            print("%4s%7d%5s%4s%2s%4d%12.3f%8.3f%8.3f%6.2f%20.3f" % ("ATOM", i+1, "C", "C", "A", i+1, sum(x_peratom[count:count+int(atom_amount[i])])/int(atom_amount[i]), sum(y_peratom[count:count+int(atom_amount[i])])/int(atom_amount[i]), sum(z_peratom[count:count+int(atom_amount[i])])/int(atom_amount[i]),1.00,(low_affinity_amount[b]+high_affinity_amount[b])/2), file = dummy)
            b += 1
            print("TER", file = dummy)
        else:
            print("%4s%7d%5s%4s%2s%4d%12.3f%8.3f%8.3f%6.2f%20.3f" % ("ATOM", i+1, "C", "C", "A", i+1, sum(x_peratom[count:count+int(atom_amount[i])])/int(atom_amount[i]), sum(y_peratom[count:count+int(atom_amount[i])])/int(atom_amount[i]), sum(z_peratom[count:count+int(atom_amount[i])])/int(atom_amount[i]),1.00,100000000000000), file = dummy)
            print("TER", file = dummy)
        count += int(atom_amount[i])







