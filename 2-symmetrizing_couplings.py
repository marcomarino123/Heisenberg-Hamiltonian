import re
import numpy
import os

###considering a supercell or not
supercell=False

###VARIABLES
file_coupling="couplings.data"
file_coordinates="atoms.coordinates.data"
file_magnetic_moments="atoms_magnetic.data"
max_distance_coupling=6.0
# here you have to indicate the number of shells you want to consider
shells={0: '3', 1: '5', 2: '6'}
# layer 0 is the molecule layer, the other layers are the slab layers
# the item in the dictionary is the minimum of the respective layer 
layers={0: '7.0', 1: '5.0', 2: '3.0', 3: '0.0'}
##the 0 element has to be the higher along z direction

###FORMAT VARIABLES
# format of the couplings.data file | Atom1 Atom2 | relative coordinates of the cell in which are the atoms rx ry rz | distance between atoms | coupling
# format of the atoms_coordinates.data | x y z

couplings_file=numpy.loadtxt(file_coupling)[:,:]
couplings_file[couplings_file[:,5].argsort()]

#saving all the couplings inside the max_distance value
count=0
for i in range(len(couplings_file[:,0])):
    if (couplings_file[i,5]<max_distance_coupling):
        count=count+1
couplings_array=numpy.zeros((count,7))
for i in range(0,count):
    couplings_array[i,:]=couplings_file[i,:]

#variable counting number atoms per layer
count_number_atoms_per_layer=numpy.zeros((len(layers),2),dtype=int)
atoms_coordinates_file=numpy.loadtxt(file_coordinates)[:,:]
atoms_coordinates_file[(-atoms_coordinates_file[:,2]).argsort()]

for key,value in layers.items():
    if key == 0:
        for i in range(len(atoms_coordinates_file[:,0])):
            if(atoms_coordinates_file[i,2]>float(value)):
                count_number_atoms_per_layer[key,1]=count_number_atoms_per_layer[key,1]+1
    else:
        for i in range(len(atoms_coordinates_file[:,0])):
            if(atoms_coordinates_file[i,2]>float(value) and atoms_coordinates_file[i,2]<=float(preceding_value)):
                count_number_atoms_per_layer[key,1]=count_number_atoms_per_layer[key,1]+1
    preceding_value=value
##writing counting number atoms as an array with: | layer i | first atom | last atom |
preceding_value=0
for key,value in layers.items():
    count_number_atoms_per_layer[key,0]=preceding_value
    count_number_atoms_per_layer[key,1]=count_number_atoms_per_layer[key,1]+preceding_value-1
    preceding_value=count_number_atoms_per_layer[key,1]+1

###TEST printing
for key,value in layers.items():
   print(key,count_number_atoms_per_layer[key,0],count_number_atoms_per_layer[key,1])

####isolating the couplings in the different layers
for key,value in layers.items():
    file=open(f"{key}.couplings.data","w+")
    for i in range(len(couplings_array[:,0])):
        if(couplings_array[i,0]>=count_number_atoms_per_layer[key,0] and couplings_array[i,0]<count_number_atoms_per_layer[key,1] and couplings_array[i,1]>=count_number_atoms_per_layer[key,0] and couplings_array[i,1]<count_number_atoms_per_layer[key,1]):
            file.write(' '.join([str(element) for element in couplings_array[i,:]]) + "\n")
    file.close()

####isolating the couplings between layers
list_layers=[key for key in layers.keys()]
def k_permutations(items, n):
    if n==0: 
        yield []
    else:
        for item in items:
            for kp in k_permutations(items, n-1):
                if item not in kp:
                    yield [item] + kp
####selecting only adjacent layers, and writing the respective couplings in the respective files, imposing order (the coplings the other way around are obtained through s_permutations
for kp in k_permutations(list_layers, 2):
    if kp[0]<kp[1]: 
        file=open(f"{kp[0]}{kp[1]}.couplings.data","w+")
        for i in range(len(couplings_array[:,0])):
            s=[0,1]
            for s_t in k_permutations(s,2):
                if(couplings_array[i,0]>=count_number_atoms_per_layer[kp[s_t[0]],0] and couplings_array[i,0]<count_number_atoms_per_layer[kp[s_t[0]],1] and couplings_array[i,1]>=count_number_atoms_per_layer[kp[s_t[1]],0] and couplings_array[i,1]<count_number_atoms_per_layer[kp[s_t[1]],1]):
                    file.write(' '.join([str(element) for element in couplings_array[i,:]]) + "\n")
        file.close()
#
##### the format of file: atom | atom | R (distance between cells)| distance atom-atom | coupling
#### in supercell calculations the coupling between nearby supercells is zero, so in order to avoid spurious effects module(R)!=0 results can be excluded
def function_mean_couplings(shells, file, supercell):
    count1=0
    file[file[:,5].argsort()]
    numb=numpy.zeros(len(shells))
    coupling_average=numpy.zeros(len(shells))
    while count1<int(len(file[:,0])):
        count2=0
        flag2=1
        while flag2!=0:
            if count2<int(len(shells)):
                if file[count1,5]<float(shells[count2]):
                    if supercell==True:
                        moduleR=numpy.sqrt(pow(file[count1,2],2)+pow(file[count1,3],2)+pow(file[count1,4],2))
                        if moduleR==0:
                            coupling_average[count2]=coupling_average[count2]+file[count1,6]
                            numb[count2]=numb[count2]+1
                            flag2=0
                            count1=count1+1
                        else:
                            count2=count2+1
                    else:
                        coupling_average[count2]=coupling_average[count2]+file[count1,6]
                        numb[count2]=numb[count2]+1
                        flag2=0
                        count1=count1+1
                else:
                    count2=count2+1
            else:
                flag2=0
                count1=count1+1
    for r in range(len(shells)):
        if(numb[r]!=0):
            print("                                     Shell:   "+str(r)+" Value:   "+str(coupling_average[r]/numb[r]))
            #+" Number Atoms:  "+str(numb[r]))
        else:
            print("                                     Shell:   "+str(r)+" Value:   "+str(0.0))
            #+" Number Atoms:  "+str(numb[r]))

print("MAX_DISTANCE: "+str(max_distance_coupling))
print("SHELLS: "+str(shells))
print("LAYERS: "+str(layers))

####average of the different couplings in the different layers
####Here I am averaging considering a circular symmetry, i.e. I am averaging the couplings with respect to the modulus of the distance
print("INTRA-LAYERS")
coupling_average=numpy.zeros((len(layers),len(shells)))
i=0
while i < len(layers):
    if(os.stat(f"{i}.couplings.data").st_size == 0):
        print(f"File {i} Empty")
    else:
        file=numpy.loadtxt(f"{i}.couplings.data")[:,:]
        print("Layer:   "+str(i))
        function_mean_couplings(shells,file,supercell)
    i=i+1

print("INTER-LAYERS")
coupling_average=numpy.zeros((len(layers),len(shells)))
i=0
while i < len(layers)-1:
    k=i+1
    while k < len(layers):
        if(os.stat(f"{i}{k}.couplings.data").st_size == 0):
            print(f"File {i}{k} Empty")
        else:
            file=numpy.loadtxt(f"{i}{k}.couplings.data")[:,:]
            print("Layer:   "+str(i)+" "+str(k))
            function_mean_couplings(shells,file,supercell)
            print("\n")
        k=k+1
    i=i+1

import sys
numpy.set_printoptions(threshold=sys.maxsize)
print("INTRA-LAYER --> DIFFERENTIATING BETWEEN UP-UP AND UP-DOWN AND VICEVERSA")
##I am substituing the number of the atom in the file_couplings with the proper magnetization from the file_magnetic_moments
####the format of the file_magnetic_moments is (atom) magnetic moment
def function_mean_couplings_magnetic_configurations(shells, file, magnetic_moments, supercell):
    file[file[:,5].argsort()]
    numb=numpy.zeros((2,len(shells)))
    coupling_average=numpy.zeros((2,len(shells)))
    count1=0
    while count1<int(len(file[:,0])):
        count2=0
        count3=0
        s=0
        ###substituing to the atoms number their magnetic moments
        while (s < int(len(magnetic_moments)) and count2+count3 < 2):
            if(int(file[count1,0])==s and count3==0):
                file[count1,0]=magnetic_moments[s]
                count3=1
            if(int(file[count1,1])==s and count2==0):
                file[count1,1]=magnetic_moments[s]
                count2=1
            s=s+1
        if ((file[count1,0]>0 and file[count1,1]>0) or (file[count1,0]<0 and file[count1,1]<0)):
            i=0
        else:
            i=1
        flag=1
        count4=0
        while flag!=0:
            if count4<int(len(shells)):
                if file[count1,5]<float(shells[count4]):
                    if supercell==True:
                        moduleR=numpy.sqrt(pow(file[count1,2],2)+pow(file[count1,3],2)+pow(file[count1,4],2))
                        if moduleR==0:
                            coupling_average[i,count4]=coupling_average[i,count4]+file[count1,6]
                            numb[i,count4]=numb[i,count4]+1
                            flag=0
                            count1=count1+1
                        else:
                            count4=count4+1
                    else:
                        coupling_average[i,count4]=coupling_average[i,count4]+file[count1,6]
                        numb[i,count4]=numb[i,count4]+1
                        flag=0
                        count1=count1+1
                else:
                    count4=count4+1
            else:
                flag=0
                count1=count1+1
    for t in range(2):
        if t==0:
            print("FM")
        else:
            print("AF")
        for r in range(len(shells)):
            if(numb[t,r]!=0):
                print("                             Shell:   "+str(r)+" Value:   "+str(coupling_average[t,r]/numb[t,r]))
                      #+" Number Atoms:  "+str(numb[t,r]))
            else:
                print("                             Shell:   "+str(r)+" Value:   "+str(0.0))
                #+" Number Atoms:  "+str(numb[t,r]))

i=0
while i < len(layers):
    if(os.stat(f"{i}.couplings.data").st_size == 0):
        print(f"File {i} Empty")
    else:
        file=numpy.loadtxt(f"{i}.couplings.data")[:,:]
        magnetic_moments=numpy.loadtxt(file_magnetic_moments)[:]
        print("Layer:   "+str(i))
        function_mean_couplings_magnetic_configurations(shells, file, magnetic_moments, supercell) 
    i=i+1


print("INTER-LAYER --> DIFFERENTIATING BETWEEN UP-UP AND UP-DOWN AND VICEVERSA")
###I am substituing the number of the atom in the file_couplings with the proper magnetization from the file_magnetic_moments
#coupling_average_f=numpy.zeros((len(layers),len(shells)))
#coupling_average_a=numpy.zeros((len(layers),len(shells)))
i=0
while i < len(layers)-1:
    k=i+1
    while k < len(layers):
        if(os.stat(f"{i}{k}.couplings.data").st_size == 0):
            print(f"File {i}{k} Empty")
        else:
            file=numpy.loadtxt(f"{i}{k}.couplings.data")[:,:]
            file_magneticmoment=numpy.loadtxt(file_magnetic_moments)[:]
            print("Layer:   "+str(i)+" "+str(k))
            function_mean_couplings_magnetic_configurations(shells, file, magnetic_moments, supercell)
            print("\n")
        k=k+1
    i=i+1
