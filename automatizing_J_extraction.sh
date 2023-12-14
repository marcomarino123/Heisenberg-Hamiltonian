#!/bin/bash
#SBATCH --nodes=2            # Number of nodes 
#SBATCH --ntasks-per-node=2    # Number of MPI ranks per node    
#SBATCH --cpus-per-task=24      # number of HW threads per task (equal to OMP_NUM_THREADS*4)
#SBATCH --mem=375300MB
#SBATCH --time 10:00:00         # Walltime, format: HH:MM:SS
#SBATCH --account=IscrB_ORGAFINT_0 
#SBATCH -p g100_usr_prod
###SBATCH --qos=g100_qos_dbg
###SBATCH --error myJob.err      # std-error file
###SBATCH --output myJob.out     # std-output file
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marco.marino1@unimi.it
#SBATCH --job-name=lessnoise.atomic

module purge
module load profile/chem-phys
module load autoload qe/6.8
module load autoload wannier90

export OMP_NUM_THREADS=4
nodes=2
tasks=2
ntasks=$((nodes*tasks))
npool=1

mkdir tmpSCF
mkdir TB2J.varyingU.MaxLocalizing
echo -ne '\ne' | source ~/TB2J_env/bin/activate
source ~/TB2J_env/bin/activate

function function_extracting () {
	tmp_1=$(awk '/the Fermi energy is/{print NR}' $1)
	tmp_2=($tmp_1)
	awk -v tmp_3=$(echo ${tmp_2[${#tmp_2[@]}-1]}) 'NR == tmp_3 {print substr($0, 28, 9)}' $1
}

#Uparameterlist=(4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9 7.0 7.2 7.4 7.6 7.8 8.0 8.2 8.4)
Uparameterlist=(6.5 6.6 6.7 6.8 6.9 7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.8 7.9 8.0)
for U in "${Uparameterlist[@]}";
do

	declare -a calculationlist=("scf" "nscf")		
	for i in ${!calculationlist[@]};
	do
		calc=${calculationlist[$i]}
		inp=NiObulk.Uparameter$U.$calc.in
		out=NiObulk.Uparameter$U.$calc.out
		if [[ "$calc" == "scf" ]]; then
                 	stat="false"
	   		echo scf		
                else
                        stat="true"
			echo nscf
                fi

cat << EOF > $inp
&control
calculation = '$calc'
restart_mode='from_scratch',
prefix='NiO_2FU_U$U',
pseudo_dir = '../pseudoGBRV',
outdir='tmp/'
/
&system
input_dft = 'vdw-df-c09'
ibrav=  0, celldm(1)=7.876378487440209
nat=  4, ntyp= 3,
ecutwfc = 70, ecutrho = 280
starting_magnetization(1)= 0.0,
starting_magnetization(2)= 0.5,
starting_magnetization(3)=-0.5,
occupations='smearing', smearing='mv', degauss=0.01,
nspin=2,
lda_plus_u=.true., lda_plus_u_kind=1, U_projection_type='atomic', Hubbard_U(2)=$U, Hubbard_U(3)=$U, !!!! U !!!!
nbnd=29
nosym=.$stat.
noinv=.$stat.
/
&electrons
diago_full_acc=.$stat.
mixing_mode = 'plain'
mixing_beta = 0.3
conv_thr =  1.0d-6
/
&ions
/
CELL_PARAMETERS alat
0.50 0.50 1.00
0.50 1.00 0.50
1.00 0.50 0.50
ATOMIC_SPECIES
O    1.  o_pbe_v1.2.uspp.F.UPF
Ni1  1.  ni_pbe_v1.4.uspp.F.UPF
Ni2  1.  ni_pbe_v1.4.uspp.F.UPF
ATOMIC_POSITIONS {crystal}
O   0.25 0.25 0.25
O   0.75 0.75 0.75
Ni1 0.0  0.0  0.0
Ni2 0.5  0.5  0.5
EOF

cat << EOF > tmp_1.txt 
K_POINTS {automatic}
8 8 8 0 0 0
EOF

cat << EOF > tmp_2.txt
K_POINTS crystal
64
  0.00000000  0.00000000  0.00000000  1.562500e-02
  0.00000000  0.00000000  0.25000000  1.562500e-02
  0.00000000  0.00000000  0.50000000  1.562500e-02
  0.00000000  0.00000000  0.75000000  1.562500e-02
  0.00000000  0.25000000  0.00000000  1.562500e-02
  0.00000000  0.25000000  0.25000000  1.562500e-02
  0.00000000  0.25000000  0.50000000  1.562500e-02
  0.00000000  0.25000000  0.75000000  1.562500e-02
  0.00000000  0.50000000  0.00000000  1.562500e-02
  0.00000000  0.50000000  0.25000000  1.562500e-02
  0.00000000  0.50000000  0.50000000  1.562500e-02
  0.00000000  0.50000000  0.75000000  1.562500e-02
  0.00000000  0.75000000  0.00000000  1.562500e-02
  0.00000000  0.75000000  0.25000000  1.562500e-02
  0.00000000  0.75000000  0.50000000  1.562500e-02
  0.00000000  0.75000000  0.75000000  1.562500e-02
  0.25000000  0.00000000  0.00000000  1.562500e-02
  0.25000000  0.00000000  0.25000000  1.562500e-02
  0.25000000  0.00000000  0.50000000  1.562500e-02
  0.25000000  0.00000000  0.75000000  1.562500e-02
  0.25000000  0.25000000  0.00000000  1.562500e-02
  0.25000000  0.25000000  0.25000000  1.562500e-02
  0.25000000  0.25000000  0.50000000  1.562500e-02
  0.25000000  0.25000000  0.75000000  1.562500e-02
  0.25000000  0.50000000  0.00000000  1.562500e-02
  0.25000000  0.50000000  0.25000000  1.562500e-02
  0.25000000  0.50000000  0.50000000  1.562500e-02
  0.25000000  0.50000000  0.75000000  1.562500e-02
  0.25000000  0.75000000  0.00000000  1.562500e-02
  0.25000000  0.75000000  0.25000000  1.562500e-02
  0.25000000  0.75000000  0.50000000  1.562500e-02
  0.25000000  0.75000000  0.75000000  1.562500e-02
  0.50000000  0.00000000  0.00000000  1.562500e-02
  0.50000000  0.00000000  0.25000000  1.562500e-02
  0.50000000  0.00000000  0.50000000  1.562500e-02
  0.50000000  0.00000000  0.75000000  1.562500e-02
  0.50000000  0.25000000  0.00000000  1.562500e-02
  0.50000000  0.25000000  0.25000000  1.562500e-02
  0.50000000  0.25000000  0.50000000  1.562500e-02
  0.50000000  0.25000000  0.75000000  1.562500e-02
  0.50000000  0.50000000  0.00000000  1.562500e-02
  0.50000000  0.50000000  0.25000000  1.562500e-02
  0.50000000  0.50000000  0.50000000  1.562500e-02
  0.50000000  0.50000000  0.75000000  1.562500e-02
  0.50000000  0.75000000  0.00000000  1.562500e-02
  0.50000000  0.75000000  0.25000000  1.562500e-02
  0.50000000  0.75000000  0.50000000  1.562500e-02
  0.50000000  0.75000000  0.75000000  1.562500e-02
  0.75000000  0.00000000  0.00000000  1.562500e-02
  0.75000000  0.00000000  0.25000000  1.562500e-02
  0.75000000  0.00000000  0.50000000  1.562500e-02
  0.75000000  0.00000000  0.75000000  1.562500e-02
  0.75000000  0.25000000  0.00000000  1.562500e-02
  0.75000000  0.25000000  0.25000000  1.562500e-02
  0.75000000  0.25000000  0.50000000  1.562500e-02
  0.75000000  0.25000000  0.75000000  1.562500e-02
  0.75000000  0.50000000  0.00000000  1.562500e-02
  0.75000000  0.50000000  0.25000000  1.562500e-02
  0.75000000  0.50000000  0.50000000  1.562500e-02
  0.75000000  0.50000000  0.75000000  1.562500e-02
  0.75000000  0.75000000  0.00000000  1.562500e-02
  0.75000000  0.75000000  0.25000000  1.562500e-02
  0.75000000  0.75000000  0.50000000  1.562500e-02
  0.75000000  0.75000000  0.75000000  1.562500e-02
EOF

			if [[ "$calc" == "scf" ]]; then
				cat tmp_1.txt >> $inp
			else
				cat tmp_2.txt >> $inp
			fi

			mpirun \
			-np $ntasks \
			pw.x \
			-npool $npool -ndiag 1 \
			-input $inp \
			2>&1 >> $out
		
			if [[ "$calc" == "scf" ]]; then
				cp -r tmp/NiO_2FU_U$U* tmpSCF/
			fi
			if [[ "$calc" == "nscf" ]]; then
			####qua sono nel caso nscf	
			###	Fermi=$(function_extracting $out)
				Fermi=$(function_extracting NiObulk.Uparameter$U.nscf.out)
				mkdir wannier90.$U
				cd wannier90.$U/
				
				declare -a spinlist=("up" "down")
        			for j in ${!spinlist[@]};
        			do
                			spin=${spinlist[$j]}
					mkdir Calc_$spin
					cd Calc_$spin/
					if [[ $U < 6.2 ]]; then
						window_froz_max=$(echo $Fermi + 2.000 | bc);
						window_froz_min=$(echo $Fermi - 12.0000 | bc);
					elif [[ $U < 6.5 ]]; then 
						window_froz_max=$(echo $Fermi + 1.000 | bc);
                                                window_froz_min=$(echo $Fermi - 12.0000 | bc);
					else
						window_froz_max=$(echo $Fermi + 0.000 | bc);
						window_froz_min=$(echo $Fermi - 12.0000 | bc)
					fi;

cat << EOF > Calc_$spin.notcomplete.win		
num_wann = 16

dis_froz_min=$window_froz_min
dis_froz_max=$window_froz_max

spin=$spin

iprint=5

Begin Projections
O: p
Ni: d
End Projections

length_unit=bohr

dis_num_iter=3000

guiding_centres=.true.

begin unit_cell_cart
bohr
3.938189244 3.938189244 7.876378487440209
3.938189244 7.876378487440209 3.938189244
7.876378487440209 3.938189244 3.938189244
end unit_cell_cart

begin atoms_frac
O 0.25 0.25 0.25
O 0.75 0.75 0.75
Ni 0.00 0.00 0.00
Ni 0.50 0.50 0.50
end atoms_frac

mp_grid = 4 4 4

begin kpoints
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.25000000
  0.00000000  0.00000000  0.50000000
  0.00000000  0.00000000  0.75000000
  0.00000000  0.25000000  0.00000000
  0.00000000  0.25000000  0.25000000
  0.00000000  0.25000000  0.50000000
  0.00000000  0.25000000  0.75000000
  0.00000000  0.50000000  0.00000000
  0.00000000  0.50000000  0.25000000
  0.00000000  0.50000000  0.50000000
  0.00000000  0.50000000  0.75000000
  0.00000000  0.75000000  0.00000000
  0.00000000  0.75000000  0.25000000
  0.00000000  0.75000000  0.50000000
  0.00000000  0.75000000  0.75000000
  0.25000000  0.00000000  0.00000000
  0.25000000  0.00000000  0.25000000
  0.25000000  0.00000000  0.50000000
  0.25000000  0.00000000  0.75000000
  0.25000000  0.25000000  0.00000000
  0.25000000  0.25000000  0.25000000
  0.25000000  0.25000000  0.50000000
  0.25000000  0.25000000  0.75000000
  0.25000000  0.50000000  0.00000000
  0.25000000  0.50000000  0.25000000
  0.25000000  0.50000000  0.50000000
  0.25000000  0.50000000  0.75000000
  0.25000000  0.75000000  0.00000000
  0.25000000  0.75000000  0.25000000
  0.25000000  0.75000000  0.50000000
  0.25000000  0.75000000  0.75000000
  0.50000000  0.00000000  0.00000000
  0.50000000  0.00000000  0.25000000
  0.50000000  0.00000000  0.50000000
  0.50000000  0.00000000  0.75000000
  0.50000000  0.25000000  0.00000000
  0.50000000  0.25000000  0.25000000
  0.50000000  0.25000000  0.50000000
  0.50000000  0.25000000  0.75000000
  0.50000000  0.50000000  0.00000000
  0.50000000  0.50000000  0.25000000
  0.50000000  0.50000000  0.50000000
  0.50000000  0.50000000  0.75000000
  0.50000000  0.75000000  0.00000000
  0.50000000  0.75000000  0.25000000
  0.50000000  0.75000000  0.50000000
  0.50000000  0.75000000  0.75000000
  0.75000000  0.00000000  0.00000000
  0.75000000  0.00000000  0.25000000
  0.75000000  0.00000000  0.50000000
  0.75000000  0.00000000  0.75000000
  0.75000000  0.25000000  0.00000000
  0.75000000  0.25000000  0.25000000
  0.75000000  0.25000000  0.50000000
  0.75000000  0.25000000  0.75000000
  0.75000000  0.50000000  0.00000000
  0.75000000  0.50000000  0.25000000
  0.75000000  0.50000000  0.50000000
  0.75000000  0.50000000  0.75000000
  0.75000000  0.75000000  0.00000000
  0.75000000  0.75000000  0.25000000
  0.75000000  0.75000000  0.50000000
  0.75000000  0.75000000  0.75000000
end kpoints
bands_plot=true
bands_num_points=20
begin kpoint_path
 K1 0.00000 0.00000 0.00000  K2  0.00000 0.00000 0.50000
 K2 0.00000 0.00000 0.50000  K3  0.18750 -0.18750 0.50000
 K3 0.18750 -0.18750 0.50000 K4  0.34375 -0.18750 0.34375
 K4 0.34375 -0.18750 0.34375 K5  0.65625 0.18750 0.65625
 K5 0.65625 0.18750 0.65625  K6  0.50000 0.50000 0.50000
end kpoint_path
bands_plot_format=gnuplot
write_hr = .true.
write_xyz = .true.
EOF
cat << EOF > tmp_3.txt
num_bands = 29
num_iter=0
EOF
cat << EOF > tmp_4.txt
num_bands = 29
num_iter=0
EOF
cat << EOF > tmp_5.txt
restart=wannierise
num_bands = 29
num_iter=0
EOF


					orderlist=(1 2 3)
					for order in "${orderlist[@]}";
					do
						cp Calc_$spin.notcomplete.win  Calc_$spin.win
						if [ "$order" -eq 1 ]; then
	0000000						cat tmp_3.txt >> Calc_$spin.win
							
							echo "preprocessing"
							mpirun \
							-np $ntasks \
							wannier90.x \
							-pp Calc_$spin
cat << EOF > Calc_$spin.pw2wan
&inputpp
outdir='../../tmp/'
prefix='NiO_2FU_U$U'
seedname='Calc_$spin'
spin_component='$spin'
write_unk=.false.
write_amn=.true.
write_mmn=.true.
/
EOF
							echo "projecting"
							mpirun \
							-np $ntasks \
							pw2wannier90.x \
							-input Calc_$spin.pw2wan \
							2>&1 >> Calc_$spin.pw2wan.out
						elif [ "$order" -eq 2 ]; then
							cat tmp_4.txt >> Calc_$spin.win
							echo "wannierizing"
							mpirun \
							-np $ntasks \
							wannier90.x \
							Calc_$spin 

						else
							cat tmp_5.txt >> Calc_$spin.win
							echo "re-centering"
							mpirun \
                                                        -np $ntasks \
                                                        wannier90.x \
                                                        Calc_$spin
		       			fi;
		       		done;
					rm  Calc_$spin.notcomplete.win
				### exiting from Calc_up/down directories
		       	cd ../;
		       	done;
				### exiting from wannier90 directory
		       	cd ../;				
			### conclusion post-nscf calculations
			fi;
		#### conclusion scf+projwfc/nscf+wannier90 calculations
		done;
		#
		#### inizio dei calcoli TB2J
		cd TB2J.varyingU.MaxLocalizing
		mkdir Jij.$U
		cd Jij.$U
		#
		declare -a spinlist=("up" "down")
                for k in ${!spinlist[@]};
                do
                        spin=${spinlist[$k]}
			cp ../../$inp ./
			cp ../../wannier90.$U/Calc_$spin/Calc_"$spin".win ./
			cp ../../wannier90.$U/Calc_$spin/Calc_"$spin"_hr.dat ./
			cp ../../wannier90.$U/Calc_$spin/Calc_"$spin"_centres.xyz ./
		done;
		
		echo "coupling"
		mpirun\
		-np 1 \
		wann2J.py --posfile $inp --efermi $Fermi --kmesh 4 4 4 --elements Ni O --prefix_up Calc_up --prefix_down Calc_down --rcut 20 --nz 500 --emax 0.1 --np 40
	
		echo "downfolding"
		python ../../TB2J_downfold.py --inpath TB2J_results --outpath TB2J_results_downfold --metals Ni --ligands O --qpath 12 12 12	
		cd TB2J_results_downfold/
		../../../extracting.Jij.slab.sh
		echo $U >> ../../couplings.data
		awk 'NR<60{print}' couplings.data >> ../../couplings.data
		cd ../../../;
## conclusion U cycles
done;
