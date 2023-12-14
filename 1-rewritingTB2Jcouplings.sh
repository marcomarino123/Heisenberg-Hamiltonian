#!/bin/bash

### Implementation for the ISOTROPIC EXCHANGE

#### 1) Extraction of the exchange couplings from the TB2J output
#### 2) Calculation of Magnonic Bands associated with the k-path specified in the file kpoints.txt

#### the format of the kpoints.txt file is:
#### special k points | points in between  
#### Gamma	10
#### K		10
#### L		0
#### It is necessary to give the spin configuration of the atoms in the primitive cell
#### the format of the spins.txt file is: 
#### s1x s1y s1z
#### s2x s2y s2z ...

file='exchange.out'
max_distance=50.0

#### extracting crystal files 
nrstarting=$(awk '/Atom /{print NR}' $file)
nrfinishing=$(awk '/Total /{print NR}' $file)
echo | awk -v nrs=$(echo $nrstarting) -v nrf=$(echo $nrfinishing) '{line[NR] = $0} 
	END { i = nrs+1
	while(i<nrf){
		print line[i]
		i = i+1
	}		
}' $file > atoms.data
awk '{printf"%s\n",$1}' atoms.data > atoms.names.data
awk '{printf"%.5f %.5f %.5f\n",$2,$3,$4}' atoms.data > atoms.coordinates.data
awk '{printf"%.5f\n",$6}' atoms.data > atoms_magnetic.data
sed -n '/Cell (Angstrom):/{n;N;N;p}' $file > bravais.lattice.data

#### extracting exchange couplings
nrstarting=$(awk '/Exchange:/{print NR}' $file)
echo | awk -v nrs=$(echo $nrstarting) '{line[NR] = $0}
        END { i = nrs + 2
        while (i<=NR){
        	print line[i]
               	i = i+1
        }
}' $file | sed 's/[)(,:]//g' | awk 'NF>2{printf"%s %s %.5f %.5f %.5f %.5f %.5f\n",$1,$2,$3,$4,$5,$10,$6}' > couplings_tmp.data

#### substituing to elements names ordering numbers
function substituing(){
python - << EOF
import re
with open('atoms.names.data','r') as file1:
    lines_file1=file1.readlines()
elements=[]
for element in lines_file1:
    elements.append(element.replace("\n",""))

file2=open('couplings_tmp.data','r')
for line_file2 in file2:
    line_file2=line_file2.split()
    for idx2 in range(0,2):
        for idx1,item1 in enumerate(elements):
            if(item1==line_file2[idx2]):
                line_file2[idx2]=idx1
    print(line_file2)
EOF
}
### format of the couplings.data file | Atom1 Atom2 | relative coordinates of the cell in which are the atoms rx ry rz | distance between atoms | coupling
substituing | sed 's/[][",'"'"']//g' > couplings.data
rm couplings_tmp.data
echo "Remind to produce the appropriate Spin File and K_point File"

