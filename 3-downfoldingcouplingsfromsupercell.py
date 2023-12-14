import re
import numpy as np
import numpy.lib as npl
import os

### instead of the name of the atoms: numbers are substituted to the names
### a file containing the unit cell 
file_unit_cell='unit_cell_file.data'
### a file containing the supercell
file_supercell='supercell_file.data'
### file containing couplings in the supercell
#atom1 | atom2 | lattice vector | coupling | effective distance
file_couplings='couplings.data'
####where to save downfolded couplings
file_downfolded_couplings='couplings_downfolded.data'
#####where to save downfolded and averaged data
file_downfolded_couplings_average='couplings_downfolded_average.data'
#####modifying the file in order to make it readable from the C++ program
file_lattice_coupling_average='file_downfolded_couplings_average_modified'
### a matrix containing the epitaxial matrix
matrix=[[6,0],[-3,3]]
### primitive vectors of the surface unit cell
a=[[2.947,0,0],[2.947,5.894,0]]

def function_distance(atom1,atom2):
    d=0
    for i in range(0,3):
        d=d+(atom1[i]-atom2[i])**2
    d=np.sqrt(d)
    return d

def function_translate(atom1,translate_1,tranlate_2,surface_primitive_vectors):
    new_coordinates=np.zeros(3)
    for i in range(0,3):
        new_coordinates[i]=atom1[i]
        new_coordinates[i]=new_coordinates[i]+translate_1*surface_primitive_vectors[0][i]+tranlate_2*surface_primitive_vectors[1][i]
    return new_coordinates

### 3 columns: xyz coordinates
unit_cell=np.loadtxt(file_unit_cell)[:,:]
number_atoms1=unit_cell.shape[0]
supercell=np.loadtxt(file_supercell)[:,:]
number_atoms2=supercell.shape[0]
### how to build the supercell starting from the unit cell
if(matrix[0][0]>=matrix[0][1]):
    primary_vector1=0
else:
    primary_vector1=1

if(matrix[1][0]>=matrix[1][1]):
    primary_vector2=0
else:
    primary_vector2=1
factor1=matrix[0][(1-primary_vector1)]/matrix[0][primary_vector1]
factor2=matrix[1][(1-primary_vector2)]/matrix[1][primary_vector2]

list_equivalents={}
for i in range(0,number_atoms1):
    list_equivalents[str(i)]=[]

print(list_equivalents)
####saving the equivalent atoms and the respective position
epsilon=0.5
translated_coordinates=np.zeros(3)
for i in range(0,matrix[0][primary_vector1]):
    for j in range(0,matrix[1][primary_vector2]):
        for atom1 in range(0,number_atoms1):
            for xyz in range(0,3):
                translated_coordinates[xyz]=unit_cell[atom1][xyz]+i*(factor1*a[(1-primary_vector1)][xyz]+a[primary_vector1][xyz])+j*(factor2*a[(1-primary_vector2)][xyz]+a[primary_vector2][xyz])
            print(translated_coordinates)
            for atom2 in range(0,number_atoms2):
                distance=function_distance(supercell[atom2,:],translated_coordinates)
                if(distance<epsilon):
                   if(len(list_equivalents[str(atom1)])==0):
                        list_equivalents[str(atom1)]=[atom2]
                   else:
                        list_equivalents[str(atom1)]=list_equivalents[str(atom1)]+[atom2]

print(list_equivalents)
####substituing in the coupling file the supercell atomic names with the unitcell atomic names (numbers)
####following list of equivalents (list_equivalents)
couplings_supercell=np.loadtxt(file_couplings)
lattice_vector=np.zeros(3)
number_lines_couplings_supercell=couplings_supercell.shape[0]
####saving the downfolded couplings inside a file.data
with open(file_downfolded_couplings, 'w+') as file:
    for i in range(0,number_lines_couplings_supercell):
        count=0
        key_transformed=[0,0]
        for key,elements in list_equivalents.items():
            for element in elements:
                for j in range(0,2):
                    ##print(element,int(couplings_supercell[i,j]))
                    if(element==int(couplings_supercell[i,j])):
                        lattice_vector=couplings_supercell[i,2:5]
                        ###in order to avoid to take spurious couplings
                        if(function_distance(lattice_vector,lattice_vector)==0):
                            key_transformed[count]=int(key)
                            count=count+1
                if(count==2):
                    break
        tmp_variable=list([int(key_transformed[0]),int(key_transformed[1]),float(couplings_supercell[i,5]),float(couplings_supercell[i,6])])
        file.write(str(tmp_variable[0])+" "+str(tmp_variable[1])+" "+str(tmp_variable[2])+" "+str(tmp_variable[3])+" "+'\n')
        ##header=npl.format.header_data_from_array_1_0(tmp_variable)
        ##npl.format.write_array_header_1_0(file,header)
        ##tmp_variable.tofile(file)
        ##file.seek(0)
        ##header['shape']=(len(tmp_variable),*header['shape'][1:])
        ##npl.format.write_array_header_1_0(file,header)

#####averaging over the couplings between the same atoms at the same distance
##### Read file contents.
couplings_unitcell=np.loadtxt(f'{file_downfolded_couplings}')
##print(couplings_unitcell)
number_lines_couplings=couplings_unitcell.shape[0]
def k_permutations(items,n):
    if n==0:
        yield []
    else:
        for item in items:
            for kp in k_permutations(items,n-1):
                if item not in kp:
                    yield [item] + kp
####for each pair of atom in the unit cell saving all the distances with respective couplings in a single list, and then averging over the duplicates in the list
list_atoms_unitcell=[i for i in range(0,number_atoms1)]
number_pairs=0
for pair in k_permutations(list_atoms_unitcell,2):
    if pair[1]>pair[0]:
        number_pairs=number_pairs+1
####adding to the list self-interactions
number_pairs=number_pairs+number_atoms1

#print(couplings_unitcell)
pairs_couplings={}
pairs_couplings_averaged={}
####grouping all the couplings in a , each entry in the matrix is associated with a pair
counting=0
for pair in k_permutations(list_atoms_unitcell,2):
    if pair[1]>=pair[0]:
        print(pair)
        pairs_couplings_tmp=[]
        for line in range(0,number_lines_couplings):
            if pair[0]==couplings_unitcell[line,0] and pair[1]==couplings_unitcell[line,1]:
                pairs_couplings_tmp.append(list([couplings_unitcell[line,2],couplings_unitcell[line,3]]))
        pairs_couplings[counting]=pairs_couplings_tmp
        pairs_couplings_averaged[counting]=[]
        counting=counting+1
###considering self-interaction
for i in range(0,number_atoms1):
    pairs_couplings_tmp=[]
    for line in range(0,number_lines_couplings):
        if i==couplings_unitcell[line,0] and i==couplings_unitcell[line,1]:
            pairs_couplings_tmp.append(list([couplings_unitcell[line,2],couplings_unitcell[line,3]]))
        pairs_couplings[counting]=pairs_couplings_tmp
    pairs_couplings_averaged[counting]=[]
    counting=counting+1

#print(pairs_couplings[0])
print("START AVERAGE")
########averaging the coupling with respect to the distance
epsilon=0.5
for key,element in pairs_couplings.items():
    if len(element)!=0:
        print("Before average",len(pairs_couplings[key]))
        for line in range(0,len(element)):
            if (line == 0):
                min_tmp=element[0][0]
                mean_tmp=0
                mean_tmp_old=0
                distance_tmp=min_tmp
                counting=0
            #print(matrix[line,0],min_tmp,epsilon)
            if element[line][0]<=min_tmp+epsilon:
                mean_tmp=mean_tmp+element[line][1]
                distance_tmp=distance_tmp+element[line][0]
                counting=counting+1
            else:
                mean_tmp=mean_tmp/counting
                distance_tmp=distance_tmp/counting
                print(distance_tmp,mean_tmp)
                if len(pairs_couplings_averaged[key])==0:
                    pairs_couplings_averaged[key]=[distance_tmp,mean_tmp]
                else:
                    pairs_couplings_averaged[key]=pairs_couplings_averaged[key]+[distance_tmp,mean_tmp]
                distance_tmp=0
                mean_tmp=0
                counting=0
                min_tmp=element[line][0]
        print("After average",len(pairs_couplings_averaged[key]))
print(pairs_couplings_averaged)
#####writing the averages in a new txt file adding a symmetrical number of reciprocal vectors pointing in each direction
########reciprocal vectors commensurate with the distance indicated
counting=0
with open(file_downfolded_couplings_average, 'w+') as file:
    for pair in k_permutations(list_atoms_unitcell,2):
        if pair[1]>pair[0]:
            elements=pairs_couplings_averaged[counting]
            for i in range(0,len(elements)):
                if i % 2 ==0:
                    file.write(str(pair[0])+" "+str(pair[1])+" "+str(elements[i])+" "+str(elements[i+1])+"\n")
                    file.write(str(pair[1])+" "+str(pair[0])+" "+str(elements[i])+" "+str(elements[i+1])+"\n")
                    print(str(pair[0])+" "+str(pair[1])+" "+str(elements[i])+" "+str(elements[i+1]))
            counting=counting+1
    for i in range(0,number_atoms1):
        elements=pairs_couplings_averaged[counting]
        for j in range(0,len(elements)):
            if j % 2 ==0:
                file.write(str(i)+" "+str(i)+" "+str(elements[j])+" "+str(elements[j+1])+"\n")
                print(str(i)+" "+str(i)+" "+str(elements[j])+" "+str(elements[j+1]))
        counting=counting+1

####number_shells_considered
##number_shells=5
##max_distance=10.000
##counting=0
##with open(file_shells_coupling_average, 'w+') as file3:
##    for pair in k_permutations(list_atoms_unitcell,2):
##        if pair[1]>pair[0]:
##            elements=pairs_couplings_averaged[counting]
##            counting2=0
##            for i in range(0,len(elements)):
##                if elements[i]<max_distance and counting2<=number_shells:
##                    if i % 2 ==0:
##                        print(elements[i],counting2)
##                        file3.write(str(pair[0])+" "+str(pair[1])+" "+str(elements[i+1])+" "+str(counting2)+"\n")
##                        file3.write(str(pair[1])+" "+str(pair[0])+" "+str(elements[i+1])+" "+str(counting2)+"\n")
##                        counting2=counting2+1
##            counting=counting+1
##    for i in range(0,number_atoms1):
##        elements=pairs_couplings_averaged[counting]
##        counting2=0
##        for j in range(0,len(elements)):
##            if elements[j]<max_distance and counting2<=number_shells:
##                if j % 2 ==0:
##                    file3.write(str(i)+" "+str(i)+" "+str(elements[j+1])+" "+str(counting2)+"\n")
##                    counting2=counting2+1
##        counting=counting+1

###TRANSFORMING THE FILE IN SOMETHING THAT CAN BE READ FROM THE MAGNONS PROGRAM
### ATOM1 ATOM2 LATTICEVECTOR DISTANCE J(SIGNLE VALUE)
counting=0
max_distance=10.000
###be carefull having done an average over the distance between atoms, to have a correspondence here, you have to be particularly large
epsilon=0.5
###to understand how many cell i have to consider to find the proper lattice vector
mod1=np.sqrt(np.dot(a[0],a[0]))
mod2=np.sqrt(np.dot(a[1],a[1]))
if mod1<mod2:
    count_max=int(max_distance/mod1) 
else:
    count_max=int(max_distance/mod2) 
print(count_max)
with open(file_lattice_coupling_average, 'w+') as file3:
    for pair in k_permutations(list_atoms_unitcell,2):
        if pair[1]>pair[0]:
            elements=pairs_couplings_averaged[counting]
            for i in range(0,len(elements)):
                if i % 2 ==0:
                    distance=elements[i]
                    for count1 in range(-count_max,count_max,1):
                        for count2 in range(-count_max,count_max,1):
                            translated_coordinates=function_translate(unit_cell[pair[1]],count1,count2,a)
                            translated_distance=function_distance(translated_coordinates,unit_cell[pair[0]])
                            if count1==count2==0:
                                print("STATUS ",pair[0],pair[1],int(i/len(elements)*100),distance,count1,count2)
                                print("translation", translated_distance,distance)
                                ##print("translated coordinates ", unit_cell[pair[1]], translated_coordinates)
                            if abs(translated_distance-distance)<=epsilon:
                                file3.write(str(pair[0])+" "+str(pair[1])+" "+str(int(count1))+" "+str(int(count2))+" "+" 0 "+str(elements[i])+" "+str(elements[i+1])+"\n")
                                file3.write(str(pair[1])+" "+str(pair[0])+" "+str(-int(count1))+" "+str(-int(count2))+" "+" 0 "+str(elements[i])+" "+str(elements[i+1])+"\n")
            counting=counting+1
    for i in range(0,number_atoms1):
        elements=pairs_couplings_averaged[counting]
        for j in range(0,len(elements)):
            if j % 2 ==0:
                distance=elements[j]
                for count1 in range(-count_max,count_max):
                    for count2 in range(-count_max,count_max):
                        translated_coordinates=function_translate(unit_cell[i],count1,count2,a)
                        translated_distance=function_distance(translated_coordinates,unit_cell[i])
                        if abs(translated_distance-distance)<epsilon:
                            file3.write(str(i)+" "+str(i)+" "+str(int(count1))+" "+str(int(count2))+" "+" 0 "+str(elements[j])+" "+str(elements[j+1])+"\n")
                            file3.write(str(i)+" "+str(i)+" "+str(-int(count1))+" "+str(-int(count2))+" "+" 0 "+str(elements[j])+" "+str(elements[j+1])+"\n")
        counting=counting+1

