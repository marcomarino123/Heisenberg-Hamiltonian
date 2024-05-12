#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdio>
#include <complex>


using namespace std;
// declaring C function/libraries in the C++ code
extern "C"
{
// wrapper of the Fortran Lapack library into C
#include <complex.h>
#include <stdio.h>
#include <omp.h>
#include <lapacke.h>
}
const double pigreco = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

///////// START GENERIC FUNCTIONS USED DURING THE EXECUTION
double function_det(double *a, double *b, double *c)
{
	double det;
	det=0;
	det = a[0] * (b[1] * c[2] - c[1] * b[2]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
	return det;
};

double *function_vector_product(double *a, double *b)
{
	double *c;
	c=new double[3];
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
	return c;
};

complex<double> function_scalar_product_2(complex<double>** a, complex<double>** b, int l, int m, int n)
{
	complex<double> al[n];
	complex<double> bm[n];
	for(int j=0;j<n;j++){
			al[j]=a[l][j];
			bm[j]=b[m][j];
	}
	complex<double> product;
	product.real(0.0);
	product.imag(0.0);
	for (int i = 0; i < n; i++)
		product = product + conj(al[i]) * (bm[i]);
	return product;
};

double function_quantum_of_phase(double phase)
{
	/// here i am assuming as main intervall the one between -\pi and \pi
	int time = (int)(phase / (2 * pigreco));
	cout << phase << " " << time << endl;
	phase = phase - time * 2 * pigreco;
	return phase;
};

double **function_skew_symmetric(double *a)
{
	double **R;
	R = new double *[3];
	for (int i = 0; i < 3; i++)
	{
		R[i] = new double[3];
		for (int j = 0; j < 3; j++)
			R[i][j] = 0.0;
	}
	R[0][1] = -a[2];
	R[0][2] = a[1];
	R[1][0] = a[2];
	R[1][2] = -a[0];
	R[2][0] = -a[2];
	R[2][1] = a[1];
	return R;
};
double function_scalar_product(double* a, double *b)
{
	double c=0;
	for (int i = 0; i < 3; i++)
		c = c + a[i] * b[i];
	return c;
};
complex<double> function_determinant(complex<double> **U, int n)
{
	/// It is convenient to diagonalize the matrix at that point the determinant is the product of the eigenvalues
	int INFO;
	char UPLO = 'L';
	char JOBVS = 'N';
	char SORT = 'N';
	char JOBZ = 'N';
	int LDA = n;
	int N = n;
	int matrix_layout = 101;
	lapack_complex_double *WR;
	int SDIM;
	int LDVS = n;
	WR = (lapack_complex_double *)malloc(n * sizeof(lapack_complex_double));
	lapack_complex_double *VS;
	VS = (lapack_complex_double *)malloc(n * n * sizeof(lapack_complex_double));
	lapack_complex_double *temporary;
	temporary = (lapack_complex_double *)malloc(n * n * sizeof(lapack_complex_double));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			temporary[i * n + j] = U[i][j].real() + _Complex_I * (U[i][j].imag());
	LAPACK_Z_SELECT1 SELECT;

	INFO = LAPACKE_zgees(matrix_layout, JOBVS, SORT, SELECT, N, temporary, LDA, &SDIM, WR, VS, LDVS);

	free(temporary);
	free(VS);

	lapack_complex_double determinant = 1.0 + _Complex_I * 0.0;

	for (int i = 0; i < n; i++)
		determinant = determinant * WR[i];

	free(WR);

	return determinant;
};

// Merges two subarrays of array[].
// First subarray is arr[begin..mid]
// Second subarray is arr[mid+1..end]
void merge(double array[], int const left,
		   int const mid, int const right)
{
	auto const subArrayOne = mid - left + 1;
	auto const subArrayTwo = right - mid;

	// Create temp arrays
	auto *leftArray = new double[subArrayOne],
		 *rightArray = new double[subArrayTwo];

	// Copy data to temp arrays leftArray[]
	// and rightArray[]
	for (auto i = 0; i < subArrayOne; i++)
		leftArray[i] = array[left + i];
	for (auto j = 0; j < subArrayTwo; j++)
		rightArray[j] = array[mid + 1 + j];

	// Initial index of first sub-array
	// Initial index of second sub-array
	auto indexOfSubArrayOne = 0,
		 indexOfSubArrayTwo = 0;

	// Initial index of merged array
	int indexOfMergedArray = left;

	// Merge the temp arrays back into
	// array[left..right]
	while (indexOfSubArrayOne < subArrayOne &&
		   indexOfSubArrayTwo < subArrayTwo)
	{
		if (leftArray[indexOfSubArrayOne] <=
			rightArray[indexOfSubArrayTwo])
		{
			array[indexOfMergedArray] =
				leftArray[indexOfSubArrayOne];
			indexOfSubArrayOne++;
		}
		else
		{
			array[indexOfMergedArray] =
				rightArray[indexOfSubArrayTwo];
			indexOfSubArrayTwo++;
		}
		indexOfMergedArray++;
	}

	// Copy the remaining elements of
	// left[], if there are any
	while (indexOfSubArrayOne < subArrayOne)
	{
		array[indexOfMergedArray] =
			leftArray[indexOfSubArrayOne];
		indexOfSubArrayOne++;
		indexOfMergedArray++;
	}

	// Copy the remaining elements of
	// right[], if there are any
	while (indexOfSubArrayTwo < subArrayTwo)
	{
		array[indexOfMergedArray] =
			rightArray[indexOfSubArrayTwo];
		indexOfSubArrayTwo++;
		indexOfMergedArray++;
	}
}

// begin is for left index and end is
// right index of the sub-array
// of arr to be sorted */
void mergeSort(double array[], int const begin, int const end)
{
	// Returns recursively
	if (begin >= end)
		return;

	auto mid = begin + (end - begin) / 2;
	mergeSort(array, begin, mid);
	mergeSort(array, mid + 1, end);
	merge(array, begin, mid, end);
}

double function_distance(double *coordinates1, double *coordinates2)
{
	double module_tmp;
	module_tmp = sqrt(pow(coordinates1[0] - coordinates2[0], 2) + pow(coordinates1[1] - coordinates2[1], 2) + pow(coordinates1[2] - coordinates2[2], 2));
	return module_tmp;
};

double function_module_matrix(double **vector)
{
	double c;
	c = sqrt(pow(vector[0][0], 2) + pow(vector[1][1], 2) + pow(vector[2][2], 2));
	return c;
};

double function_module(double *vector)
{
	double c;
	c = sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
	return c;
};

double *function_translate(double *coordinates1, double *count, double **bravais_lattice)
{
	double *coordinates_translated;
	coordinates_translated = new double[3];
	for (int i = 0; i < 3; i++)
	{
		coordinates_translated[i] = coordinates1[i];
		for (int j = 0; j < 3; j++)
			coordinates_translated[i] = coordinates_translated[i] + count[j] * bravais_lattice[j][i];
	}
	return coordinates_translated;
};
///////// END GENERIC FUNCTIONS USED DURING THE EXECUTION


///////// START DATA_STRUCTURE USED IN THE CREATION OF TABLES AND LISTS
class Node_list
{
private:
	int site;
	double **value;
	double *distance;
	Node_list *next;
public:
	Node_list()
	{
		site = -1;
		value = nullptr;
		distance = nullptr;
		next = nullptr;
	};
	void fill_node(double **value_tmp, double *distance_tmp, int site_tmp)
	{
		site = site_tmp;
		value = value_tmp;
		distance = distance_tmp;
		next=nullptr;
	};
	Node_list* push_value(double **value_tmp, double *distance_tmp)
	{
		Node_list* next_tmp=new Node_list();
		next_tmp->fill_node(value_tmp,distance_tmp,site+1);
		next=next_tmp;
		return next;
	};
	double **pull_value()
	{
		return value;
	};
	int pull_site()
	{
		return site;
	};
	double *pull_distance()
	{
		return distance;
	};
	Node_list *pull_next()
	{
		return next;
	};
	~Node_list(){
		Node_list* next_tmp;
		Node_list* prec_tmp;
		prec_tmp=this;
		next_tmp=next;
		while(next_tmp!=NULL){
			delete[] prec_tmp->pull_value();
			delete[] prec_tmp->pull_distance();
			prec_tmp=prec_tmp->pull_next();
			next_tmp=prec_tmp->pull_next();
		}
		delete[] prec_tmp;
	};
};
///////// START DATA_STRUCTURE USED IN THE CREATION OF TABLES AND LISTS

///////// START DEFINITION CLASS SPIN-OPERATOR
class spin_operator
{
private:
	/// flag to distinguish between magnetic and no-magnetic atoms
	int flag;
	/// site of the spin
	int site;
	/// coordinate of the spin
	double *site_coordinates;
	/// modulus of the spin vector
	double modulus;
	/// local frame
	double *orientation;
	complex<double> u[3];

public:
	spin_operator()
	{
		site_coordinates = nullptr;
		flag = 0;
		site = 0;
		modulus = 0;
		orientation = nullptr;
	}
	void push_values(int site_tmp, double *site_coordinates_tmp, double *spin_value_tmp);
	void define_local_frame();
	void print();
	double pull_modulus()
	{
		return modulus;
	}
	complex<double> *pull_u()
	{
		return u;
	}
	double *pull_orientation()
	{
		return orientation;
	}
	~spin_operator(){
		delete[] site_coordinates;
		delete[] orientation;
	}
};

void spin_operator::push_values(int site_tmp, double *site_coordinates_tmp, double *spin_value_tmp)
{
	site = site_tmp;
	site_coordinates = site_coordinates_tmp;
	orientation = spin_value_tmp;
	for (int i = 0; i < 3; i++)
	{
		modulus = modulus + pow(spin_value_tmp[i], 2);
	}
	modulus = sqrt(modulus);
	if (modulus != 0)
	{
		for (int j = 0; j < 3; j++)
			orientation[j] = orientation[j] / modulus;
		flag = 1;
	}
	define_local_frame();
};

void spin_operator::define_local_frame()
{
	/// the local frame can be defined only for a magnetic atom
	if (modulus != 0)
	{
		double norm1;
		double norm2;
		int counting;
		double* Real_u;
		double* Imag_u;
		/// defining lab frame
		double* z=new double[3];
		z[0] = 0;
		z[1] = 0;
		z[2] = 1;
		norm1 = function_scalar_product(z, orientation);
		/// if z lab frame is not parallel to spin orientation
		if (abs(norm1) != 1)
		{
			double* v=new double[3];
			v[0] = -orientation[1] * z[2] + orientation[2] * z[1];
			v[1] = orientation[0] * z[2] - orientation[2] * z[0];
			v[2] = -orientation[0] * z[1] + orientation[1] * z[0];
			norm2 = 0;
			for (int i = 0; i < 3; i++)
			{
				norm2 = norm2 + orientation[i] * z[i];
			}
			double** R=new double*[3];
			double** A2=new double*[3];
			for(int i=0;i<3;i++){
				R[i]=new double[3];
				A2[i]=new double[3];
			}
			double **A = function_skew_symmetric(v);
			
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					A2[i][j] = 0;
					for (int k = 0; k < 3; k++)
						A2[i][j] = A2[i][j] + A[i][k] * A[k][j];
				}
			}
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					if (i == j)
						R[i][j] = 1 + A2[i][j] / (1 + norm2);
					else
						R[i][j] = A[i][j] + A2[i][j] / (1 + norm2);
				}
			for (int i = 0; i < 3; i++){
				delete[] (double*)A[i];
				delete[] (double*)R[i];
				delete[] (double*)A2[i];
			}
			delete[] (double**)A;
			delete[] (double**)R;
			delete[] (double**)A2;
			delete[] (double*)v;
			/// logical condition already satisfied: in order to avoid the definition of a new variable let use norm1
			norm1 = 0;
			norm2 = 0;
			/// normalizing the columns in order to define the versors orthogonal to the spin orientation
			for (int i = 0; i < 3; i++)
			{
				norm1 = norm1 + pow(R[i][0], 2);
				norm2 = norm2 + pow(R[i][1], 2);
			}
			norm1 = sqrt(norm1);
			norm2 = sqrt(norm2);
			for (int i = 0; i < 3; i++)
			{
				u[i].real(R[i][0] / norm1);
				u[i].imag(R[i][1] / norm2);
			}
		}
		else
		{
			/// in the case z and orientation are parallel
			Real_u=new double[3];
			Imag_u=new double[3];
			for (int i = 0; i < 3; i++)
			{
				Real_u[i] = 0.0;
				Imag_u[i] = 0.0;
			}
			Real_u[0] = 1.0;
			if (orientation[2] < 0)
				Imag_u[1] = -1.0;
			else
				Imag_u[1] = 1.0;
			for (int i = 0; i < 3; i++)
			{
				u[i].real(Real_u[i]);
				u[i].imag(Imag_u[i]);
			}
		}
		delete[] (double*)z;
		delete[] (double*)Real_u;
		delete[] (double*)Imag_u;
	}
	else
	{
		for (int j = 0; j < 3; j++)
		{
			u[j].real(0.0);
			u[j].imag(0.0);
		}
	}
};

void spin_operator::print()
{
	cout << "Site: " << site << endl;
	for (int i = 0; i < 3; i++)
		cout << site_coordinates[i] << " ";
	cout << endl;
	cout << "Modulus: " << modulus << endl;
	cout << "Orientation:" << endl;
	for (int j = 0; j < 3; j++)
		cout << orientation[j] << " ";
	cout << endl;
	cout << "Local frame:" << endl;
	for (int j = 0; j < 3; j++)
		cout << u[j] << endl;
};
///////// END DEFINITION CLASS SPIN-OPERATOR

///////// START DEFINITION CLASS LATTICE OF THE SYSTEM
class lattice_crystal
{
private:
	int number_atoms;
	double **atoms_coordinates;
	double **bravais_lattice;
	double volume;

public:
	lattice_crystal()
	{
		number_atoms = 0;
		atoms_coordinates = nullptr;
		bravais_lattice = nullptr;
		volume = 0.0;
	};
	void push_values(ifstream *bravais_lattice_file, ifstream *atoms_coordinates_file);
	void print();
	int pull_number_atoms()
	{
		return number_atoms;
	}
	double *pull_sitei_coordinates(int site_tmp)
	{
		return atoms_coordinates[site_tmp];
	}
	double **pull_bravais_lattice()
	{
		return bravais_lattice;
	}
	double **pull_atoms_coordinates()
	{
		return atoms_coordinates;
	}
	double pull_volume()
	{
		return volume;
	}
};

void lattice_crystal::push_values(ifstream *bravais_lattice_file, ifstream *atoms_coordinates_file)
{
	double vec1[3];
	double vec2[3];
	double vec3[3];
	bravais_lattice_file->seekg(0);
	atoms_coordinates_file->seekg(0);
	bravais_lattice = new double *[3];
	for (int i = 0; i < 3; i++)
	{
		bravais_lattice[i] = new double[3];
		for (int j = 0; j < 3; j++)
			*bravais_lattice_file >> bravais_lattice[i][j];
	}
	string line_tmp;
	while (atoms_coordinates_file->peek() != EOF)
	{
		getline(*atoms_coordinates_file, line_tmp);
		number_atoms++;
	}
	atoms_coordinates = new double *[number_atoms];
	atoms_coordinates_file->seekg(0);
	for (int i = 0; i < number_atoms; i++)
	{
		atoms_coordinates[i] = new double[3];
		for (int j = 0; j < 3; j++)
			*atoms_coordinates_file >> atoms_coordinates[i][j];
	}
	for (int i = 0; i < 3; i++)
	{
		vec1[i] = bravais_lattice[0][i];
		vec2[i] = bravais_lattice[1][i];
		vec3[i] = bravais_lattice[2][i];
	}
	volume = function_det(vec1, vec2, vec3);
};

void lattice_crystal::print()
{
	cout << "Bravais Lattice:" << endl;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			cout << bravais_lattice[i][j] << " ";
		cout << endl;
	}
	cout << "Atoms coordinates:" << endl;
	for (int i = 0; i < number_atoms; i++)
	{
		for (int j = 0; j < 3; j++)
			cout << atoms_coordinates[i][j] << " ";
		cout << endl;
	}
};
///////// END DEFINITION CLASS LATTICE OF THE SYSTEM

///////// START DEFINITION CLASS SPIN-LATTICE OF THE SYSTEM
class lattice_spin
{
private:
	int number_atoms;
	lattice_crystal *lattice;
	spin_operator *spin;
	int number_magnetic_atoms;
	double **bravais_lattice;
	double **magnetic_atom_coordinates;
	int *ordering_from_magnetic_to_crystal;
	int *ordering_from_crystal_to_magnetic;
public:
	lattice_spin()
	{
		number_atoms = 0;
		lattice = nullptr;
		spin = nullptr;
		number_magnetic_atoms = 0;
		magnetic_atom_coordinates=nullptr;
		bravais_lattice=nullptr;
		ordering_from_magnetic_to_crystal=nullptr;
		ordering_from_crystal_to_magnetic=nullptr;
	}
	void push_values(lattice_crystal *lattice_tmp, ifstream *spin_atoms_file);
	void print();
	spin_operator *pull_spin()
	{
		return spin;
	}
	int pull_number_atoms()
	{
		return number_atoms;
	}
	int pull_number_magnetic_atoms()
	{
		return number_magnetic_atoms;
	}
	double** pull_bravais_lattice()
	{
		return bravais_lattice;
	}
	double **pull_magnetic_atoms_coordinates()
	{
		return magnetic_atom_coordinates;
	}
	int *pull_ordering_from_magnetic_to_crystal()
	{
		return ordering_from_magnetic_to_crystal;
	}
	int *pull_ordering_from_crystal_to_magnetic()
	{
		return ordering_from_crystal_to_magnetic;
	}
};

void lattice_spin::push_values(lattice_crystal *lattice_tmp, ifstream *spin_atoms_file)
{
	lattice = lattice_tmp;
	number_atoms = lattice->pull_number_atoms();
	bravais_lattice = lattice->pull_bravais_lattice();
	
	int ordering_tmp[number_atoms];
	double spin_value_tmp[number_atoms][3];
	ordering_from_crystal_to_magnetic=new int[number_atoms];
	number_magnetic_atoms=0;

	spin_atoms_file->seekg(0);
	for (int i= 0;i<number_atoms;i++)
	{
		for (int j=0;j<3;j++){
			*spin_atoms_file >> spin_value_tmp[i][j];
		}

		if (function_module(spin_value_tmp[i])!= 0){
			ordering_tmp[number_magnetic_atoms]=i;
			ordering_from_crystal_to_magnetic[i]=number_magnetic_atoms;
			number_magnetic_atoms = number_magnetic_atoms + 1;
		}
		else
			ordering_from_crystal_to_magnetic[i]=number_atoms;
	}


	magnetic_atom_coordinates=new double*[number_magnetic_atoms];
	ordering_from_magnetic_to_crystal=new int[number_magnetic_atoms];
		for (int i=0;i<number_magnetic_atoms;i++){
		ordering_from_magnetic_to_crystal[i]=ordering_tmp[i];
	}
	
	spin=new spin_operator[number_magnetic_atoms];

	double* spin_value;

	for (int i=0;i<number_magnetic_atoms;i++)
	{
		magnetic_atom_coordinates[i]=new double[3];
		magnetic_atom_coordinates[i]=lattice->pull_sitei_coordinates(ordering_from_magnetic_to_crystal[i]);
		spin_value=new double[3];
		for (int r = 0; r<3; r++)
			spin_value[r]=spin_value_tmp[ordering_from_magnetic_to_crystal[i]][r];
		spin[i].push_values(ordering_from_magnetic_to_crystal[i],lattice->pull_sitei_coordinates(ordering_from_magnetic_to_crystal[i]),spin_value);
	}
	cout<<" from magnetic to crystal"<<endl;
	for (int i=0;i<number_magnetic_atoms;i++)
		cout<<i<<" "<<ordering_from_magnetic_to_crystal[i]<<endl;
	cout<<" from crystal to magnetic"<<endl;
	for (int i=0;i<number_atoms;i++)
		cout<<i<<" "<<ordering_from_crystal_to_magnetic[i]<<endl;
};

void lattice_spin::print()
{
	lattice->print();
	for (int i=0;i<number_magnetic_atoms;i++)
		spin[i].print();
};
///////// END DEFINITION CLASS SPIN-LATTICE OF THE SYSTEM


///////// START DEFINITION CLASS J-LATTICE OF THE SYSTEM
class lattice_J
{
private:
	/// maximum distance considered in the couplings between atoms
	double max_distance_coupling;
	/// two atoms i and j can be coupled at different distances
	/// to save the different coupling with respect to the relative distance
	/// between the atoms i and j a dynamical list is used
	Node_list ***Jij_lattice;
	lattice_spin *spin_lattice;

public:
	lattice_J()
	{
		max_distance_coupling = 0;
		Jij_lattice = nullptr;
		spin_lattice = nullptr;
	}
	void push_values(int max_distance_coupling_tmp, lattice_spin *spin_lattice_tmp, ifstream *couplings_file);
	void building_values(int max_distance_coupling_tmp, lattice_spin *spin_lattice_tmp);
	void print();
	/// extracting the FT of the coupling ij
	complex<double> ****pull_fft(double *kpoint_tmp);
	void print_fft(double *kpoint_tmp);
	lattice_spin *pull_spin_lattice()
	{
		return spin_lattice;
	}
	void free_fft(complex<double> ****);
	///// here the system has no spirals, so the matrix rotation multiplied to the Jij is the identity
};

void lattice_J::push_values(int max_distance_coupling_tmp, lattice_spin *spin_lattice_tmp, ifstream *couplings_file)
{
	spin_lattice = spin_lattice_tmp;
	int number_magnetic_atoms = spin_lattice->pull_number_magnetic_atoms();
	max_distance_coupling = max_distance_coupling_tmp;
	int number_atoms = spin_lattice->pull_number_atoms();

	double distance_tmp[3];
	double *distance;
	double distance_module;
	double **coupling;
	double **coupling_tmp_allocate1;
	double *coupling_tmp_allocate2;
	double coupling_tmp;
	double** temp;

	/// initialization heads of the lattice_J
	Jij_lattice=new Node_list**[number_magnetic_atoms];
	for (int i=0;i<number_magnetic_atoms;i++){
		Jij_lattice[i]=new Node_list*[number_magnetic_atoms];
		for (int j=0;j<number_magnetic_atoms;j++)
			Jij_lattice[i][j]=new Node_list();
	}
	/// pointers to go through the dynamical lists for each ij
	Node_list ***Jij_lattice_next;
	Jij_lattice_next = new Node_list**[number_magnetic_atoms];
	for (int i=0;i<number_magnetic_atoms;i++)
		Jij_lattice_next[i] = new Node_list*[number_magnetic_atoms];

	couplings_file->seekg(0);
	int i; int j; int site;
	int i_tmp; int j_tmp;
	while (couplings_file->peek()!=EOF)
	{
		distance_module = 0.0;
		*couplings_file>>i_tmp;
		*couplings_file>>j_tmp;
		/// distance between the cells of the two atoms in relative coordinates
		for (int l=0;l<3;l++)
			*couplings_file >> distance_tmp[l];
		*couplings_file >> distance_module >> coupling_tmp;
		if (distance_module <= max_distance_coupling){
			///convertion from crystal ordering of the atoms, to spin_crystal ordering
			i=spin_lattice->pull_ordering_from_crystal_to_magnetic()[i_tmp];
			j=spin_lattice->pull_ordering_from_crystal_to_magnetic()[j_tmp];
			///allocating_memory
			distance=new double[3];
			for(int s=0;s<3;s++)
				distance[s]=distance_tmp[s];
			/// initialization values to push into Jij (for each distance)
			/// considering only the isotrpic part of the coupling and putting it on the diagonal of a 3x3 matrix
			/// the generalization to anisotropic terms is quite straightforward
			coupling=new double*[3];
			for(int s=0;s<3;s++)
				coupling[s]=new double[3];
			for(int q=0;q<3;q++){
				if(q==0){
					coupling[0][0]=-coupling_tmp;
					cout<<coupling[0][0]<<endl;
				}
				for(int t=0;t<3;t++){
					if(t==q)
						coupling[q][t]=coupling[0][0];
					else
						coupling[q][t]=0.0;
				}
			}
			site = Jij_lattice[i][j]->pull_site();
			if(site < 0)
			{
				Jij_lattice[i][j]->fill_node(coupling,distance,0);
				Jij_lattice_next[i][j]=Jij_lattice[i][j];
			}
			else
				Jij_lattice_next[i][j]=Jij_lattice_next[i][j]->push_value(coupling, distance);
		}
		cout<<i<<" "<<j<<" "<<site<<" "<<distance_module<<endl;
		cout<<Jij_lattice_next[i][j]<<endl;
	}

	cout<<"checking now atoms with no entries"<<endl;
	///// initializing the atoms not having any coupling to coupling 0
	for (int i = 0; i < number_magnetic_atoms; i++)
		for (int j = 0; j < number_magnetic_atoms; j++)
		{
			cout<<i<<" "<<j<<endl;
			if(Jij_lattice[i][j]!=NULL){
				site = Jij_lattice[i][j]->pull_site();
				cout<<site<<endl;
				if(site < 0)
				{
					coupling=new double*[3];
					distance=new double[3];
					for (int r=0;r<3;r++)
					{
						distance[r] = 0.0;
						coupling[r]=new double[3];
						for (int t=0;t<3;t++)
							coupling[r][t] = 0.0;
					}
					Jij_lattice[i][j]->fill_node(coupling,distance,0);
				}
				cout<<"next pair "<<endl;
			}
		}
		
	for(int i=0;i<number_magnetic_atoms;i++)
		delete[] (double*) Jij_lattice_next[i];
	delete[] (double**) Jij_lattice_next;
	
};

void lattice_J::building_values(int max_distance_coupling_tmp, lattice_spin *spin_lattice_tmp)
{
	spin_lattice=spin_lattice_tmp;
	double **bravais_lattice = spin_lattice->pull_bravais_lattice();

	/// understanding how many cell to consider in the supercell building
	double min_module = 100;
	double module_tmp = 0;
	for (int i = 0; i < 3; i++)
	{
		module_tmp = sqrt(pow(bravais_lattice[i][0], 2) + pow(bravais_lattice[i][1], 2) + pow(bravais_lattice[i][2], 2));
		if (module_tmp < min_module)
			min_module = module_tmp;
	}
	int count_max;
	/// maybe truncation should be done for excess
	count_max = int(max_distance_coupling_tmp / min_module) + 1;
	//// we can build the supercell filtering the atoms which are too far from the primitive cell atoms
	int number_magnetic_atoms = spin_lattice->pull_number_magnetic_atoms();
	double **magnetic_atoms_coordinates = spin_lattice->pull_magnetic_atoms_coordinates();

	/// initialization heads of the lattice_J
	Jij_lattice = new Node_list **[number_magnetic_atoms];
	for (int i=0;i<number_magnetic_atoms;i++)
	{
		Jij_lattice[i] = new Node_list *[number_magnetic_atoms];
		for (int j=0;j<number_magnetic_atoms;j++)
			Jij_lattice[i][j] = new Node_list();
	}
	/// pointers to go through the dynamical lists for each ij
	Node_list ***Jij_lattice_next;
	Jij_lattice_next = new Node_list **[number_magnetic_atoms];
	for (int i = 0; i < number_magnetic_atoms; i++){
		Jij_lattice_next[i] = new Node_list *[number_magnetic_atoms];
		for (int j=0;j<number_magnetic_atoms;j++)
			Jij_lattice_next[i][j]=Jij_lattice[i][j];
	}

	double distance;
	int count1;
	int count2;
	int count3;
	double **couplingmatrix_tmp = new double *[3];
	for (int i=0;i<3;i++)
		couplingmatrix_tmp[i] = new double[3];

	////reading the effective couplings from the file couplings_effective.data
	/// Acthung: Symmetrized data!!!!
	//// number_shells
	//// distances of the respective shells
	/// i j J shell
	ifstream couplings_effective_file;
	couplings_effective_file.open("couplings_effective.data");
	int l_tmp;
	int m_tmp;
	int l;
	int m;
	double distance_tmp;
	double coupling_tmp;
	int number_shells;
	int shell_tmp;
	couplings_effective_file >> number_shells;
	double shells_distances_tmp[number_shells];
	for (int r = 0; r < number_shells; r++)
		couplings_effective_file >> shells_distances_tmp[r];

	double matrix_coupling_tmp[number_magnetic_atoms][number_magnetic_atoms][number_magnetic_atoms];
	for (int i=0;i<number_magnetic_atoms;i++)
		for (int j=0;j<number_magnetic_atoms;j++)
			for (int r=0; r<number_shells;r++)
				matrix_coupling_tmp[i][j][r] = 0;

	while (couplings_effective_file.peek() != EOF)
	{
		couplings_effective_file >> l_tmp >> m_tmp;
		l=spin_lattice->pull_ordering_from_crystal_to_magnetic()[l_tmp];
		m=spin_lattice->pull_ordering_from_crystal_to_magnetic()[m_tmp];
		couplings_effective_file >> coupling_tmp;
		couplings_effective_file >> shell_tmp;
		matrix_coupling_tmp[l][m][shell_tmp] = coupling_tmp;
	}
	couplings_effective_file.close();

	int flag2;
	int count_shell;
	int flag1;
	int site;
	int flag;
	double *magnetic_coordinates_translated;
	double *count = new double[3];
	double *count4;
	for (int i = 0; i < number_magnetic_atoms; i++)
	{
		for (int j = 0; j < number_magnetic_atoms; j++)
		{
			for (int count1 = -count_max; count1 <= count_max; count1++)
				for (int count2 = -count_max; count2 <= count_max; count2++)
					for (int count3 = -count_max; count3 <= count_max; count3++)
						/// avoiding self-interaction
						if ((count1 != 0) || (count2 != 0) || (count3 != 0))
						{
							count[0] = double(count1);
							count[1] = double(count2);
							count[2] = double(count3);
							magnetic_coordinates_translated = function_translate(magnetic_atoms_coordinates[j], count, bravais_lattice);
							distance = function_distance(magnetic_atoms_coordinates[i], magnetic_coordinates_translated);
							delete[] magnetic_coordinates_translated;
							cout << i << " " << j << " " << count[0] << " " << count[1] << " " << count[2] << " " << distance << " " << magnetic_coordinates_translated[0] << " " << magnetic_coordinates_translated[1] << " " << magnetic_coordinates_translated[2] << endl;
							cout << "i " << magnetic_atoms_coordinates[i][0] << " " << magnetic_atoms_coordinates[i][1] << " " << magnetic_atoms_coordinates[i][2] << endl;
							cout << "j " << magnetic_atoms_coordinates[j][0] << " " << magnetic_atoms_coordinates[j][1] << " " << magnetic_atoms_coordinates[j][2] << endl;
							cout << "T " << magnetic_coordinates_translated[0] << " " << magnetic_coordinates_translated[1] << " " << magnetic_coordinates_translated[2] << endl;
							if (distance < max_distance_coupling_tmp)
							{
								double *count4 = new double[3];
								count4[0] = count[0];
								count4[1] = count[1];
								count4[2] = count[2];
								cout << "ENTERING" << endl;
								count_shell = 0;
								flag2 = 0;
								flag1 = 0;
								while (flag2 != 1)
								{
									cout << flag1 << " " << flag2 << " " << count_shell << "/" << number_shells << " " << distance << " " << shells_distances_tmp[count_shell] << " " << endl;
									if (count_shell > number_shells - 1)
									{
										flag2 = 1;
										flag1 = 1;
										cout << "OVERGOING" << endl;
									}
									else if (shells_distances_tmp[count_shell] < distance)
										count_shell = count_shell + 1;
									else
										flag2 = 1;
								}
								cout << "EXITING" << endl;
								/// considering only the isotrpic part of the coupling and putting it on the diagonal of a 3x3 matrix
								/// the generalization to anisotropic terms is quite straightforward
								couplingmatrix_tmp = new double *[3];
								for (int q = 0; q < 3; q++)
								{
									couplingmatrix_tmp[q] = new double[3];
								}
								for (int l = 0; l < 3; l++)
									for (int s = 0; s < 3; s++)
									{
										if (l == s)
										{
											if (l == 0)
											{
												if (flag1 != 1)
												{
													couplingmatrix_tmp[l][s] = matrix_coupling_tmp[i][j][count_shell];
												}
												else
												{
													couplingmatrix_tmp[l][s] = 0.0;
												}
												couplingmatrix_tmp[l][s] = -couplingmatrix_tmp[l][s];
											}
											else
											{
												couplingmatrix_tmp[l][s] = couplingmatrix_tmp[0][0];
											}
										}
										else
											couplingmatrix_tmp[l][s] = 0.0;
									}
								cout << "STARTING WRITING VALUES" << endl;
								cout << "shell:" << count_shell << endl;
								cout << "count:" << count4[0] << " " << count4[1] << " " << count4[2] << " " << endl;
								site = Jij_lattice[i][j]->pull_site();
								if(site < 0)
								{
									Jij_lattice[i][j]->fill_node(couplingmatrix_tmp,count4,0);
									Jij_lattice_next[i][j]=Jij_lattice[i][j];
								}
								else
									Jij_lattice_next[i][j]=Jij_lattice_next[i][j]->push_value(couplingmatrix_tmp,count4);
							}
						}
			cout << "END OF PAIR EVALUATION (" << i << " " << j << ")" << endl;
		}
	}
	cout << "END BUILDING" << endl;
	/// initializing the atoms not having any coupling
	/// ghost coupling... in order to avoid a Node_list-matrix not equally initialized in all its entries
	for (int i = 0; i < number_magnetic_atoms; i++)
		for (int j = 0; j < number_magnetic_atoms; j++)
		{
			if(Jij_lattice_next[i][j]->pull_site()<0)
			{
				cout << "ECCOMI " << i << j << endl;
				couplingmatrix_tmp = new double *[3];
				count4 = new double[3];
				for (int i = 0; i < 3; i++)
				{
					count4[i] = 0.0;
					couplingmatrix_tmp[i] = new double[3];
					for (int j = 0; j < 3; j++)
						couplingmatrix_tmp[i][j] = 0.0;
				}
				Jij_lattice_next[i][j]->fill_node(couplingmatrix_tmp, count4,0);
			}
		}
	cout << "END THE INITIALIZATION" << endl;
};

void lattice_J::print()
{
	//cout << "ECCOMI NEL PRINTING" << endl;
	int number_magnetic_atoms = spin_lattice->pull_number_magnetic_atoms();
	Node_list ***Jij_lattice_next;
	Jij_lattice_next = new Node_list **[number_magnetic_atoms];
	for (int i = 0; i < number_magnetic_atoms; i++)
		Jij_lattice_next[i] = new Node_list *[number_magnetic_atoms];
	int count;
	//cout << "EVALUATION PAIRS" << endl;
	for (int i=0;i<number_magnetic_atoms;i++)
		for (int j=0;j<number_magnetic_atoms;j++)
		{
			count = 0;
			cout << i << " " << j << endl;
			Jij_lattice_next[i][j] = Jij_lattice[i][j];
			while (Jij_lattice_next[i][j] != nullptr)
			{
				//cout << i << j << " Distance: ";
				for (int r = 0; r < 3; r++)
					cout << Jij_lattice_next[i][j]->pull_distance()[r] << " ";
				//cout << "   Coupling: ";
				for (int r = 0; r < 3; r++)
				{
					for (int l = 0; l < 3; l++)
						if (r == l)
							cout << Jij_lattice_next[i][j]->pull_value()[r][l] << " ";
				}
				cout << endl;
				count++;
				Jij_lattice_next[i][j] = Jij_lattice_next[i][j]->pull_next();
			}
			cout << "Conteggi lettura: " << count << endl;
		}
};

complex<double> ****lattice_J ::pull_fft(double *kpoint_tmp)
{
	int number_magnetic_atoms = spin_lattice->pull_number_magnetic_atoms();
	complex<double> Jij_tmp[3][3];
	// cout<<"K point "<< endl;
	// for (int i = 0; i < 3; i++)
	//	cout<<kpoint_tmp[i]<<" ";

	complex<double> ****Jij_k;
	Jij_k = new complex<double> ***[number_magnetic_atoms];
	for (int i = 0; i < number_magnetic_atoms; i++)
	{
		Jij_k[i] = new complex<double> **[number_magnetic_atoms];
		for (int j = 0; j < number_magnetic_atoms; j++)
		{
			Jij_k[i][j] = new complex<double>*[3];
			for (int l = 0; l < 3; l++)
				Jij_k[i][j][l] = new complex<double>[3];
		}
	}

	double sum_tmp1;
	double sum_tmp2_real;
	double sum_tmp2_imag;

	typedef std::complex<double> C;

	Node_list ***Jij_lattice_next;
	Jij_lattice_next = new Node_list **[number_magnetic_atoms];
	for (int i = 0; i < number_magnetic_atoms; i++)
		Jij_lattice_next[i] = new Node_list *[number_magnetic_atoms];

	// int count;
	for (int i = 0; i < number_magnetic_atoms; i++)
		for (int j = 0; j < number_magnetic_atoms; j++)
		{
			// count=0;
			for (int l = 0; l < 3; l++)
				for (int m = 0; m < 3; m++)
				{
					Jij_tmp[l][m] = C(0.0, 0.0);
				}
			// cout<<i<<j<<endl;
			Jij_lattice_next[i][j] = Jij_lattice[i][j];
			while (Jij_lattice_next[i][j] != nullptr)
			{
				sum_tmp1 = 0;
				sum_tmp2_real = 0.0;
				sum_tmp2_imag = 0.0;
				// cout<<i<<j<<" Distance:  ";
				for (int r = 0; r < 3; r++)
				{
					sum_tmp1 = sum_tmp1 + kpoint_tmp[r] * Jij_lattice_next[i][j]->pull_distance()[r];
					// cout<<" "<<Jij_lattice_next[i][j]->pull_distance()[r];
					// cout<<"kr "<<i<<" "<<j<<" "<<sum_tmp1<<" "<<kpoint_tmp[r]<<" "<<Jij_lattice_next[i][j]->pull_distance()[r]<<endl;
				}
				sum_tmp2_real = cos(2 * pigreco * sum_tmp1);
				// cout<<"prova "<<cos(2*pigreco)<<endl;
				sum_tmp2_imag = sin(2 * pigreco * sum_tmp1);
				// cout<<"ORA  "<<Jij_lattice_next[i][j]->pull_distance()[0]<<Jij_lattice_next[i][j]->pull_distance()[1]<<Jij_lattice_next[i][j]->pull_distance()[2]<<endl;
				// cout<<"ECCO  "<<sum_tmp1<<endl;
				// cout<<" "<<sum_tmp2_real<<" +i("<<sum_tmp2_imag<<") ";
				// cout<<" Coupling:  ";
				// VALORI MANEGGIATI GLI STESSI IN INPUT: ERRORE E' NELLA TRASFORMATA DI FOURIER????
				// cout<<"Jcoupling "<<count<<" "<<Jij_lattice_next[i][j]->pull_value()[0][0]<<" "<<sum_tmp1<<endl;
				for (int l = 0; l < 3; l++)
				{
					for (int m = 0; m < 3; m++)
					{
						Jij_tmp[l][m].real(Jij_tmp[l][m].real() + Jij_lattice_next[i][j]->pull_value()[l][m] * sum_tmp2_real);
						Jij_tmp[l][m].imag(Jij_tmp[l][m].imag() - Jij_lattice_next[i][j]->pull_value()[l][m] * sum_tmp2_imag);
						// cout<<Jij_tmp[l][m]<<" ";
						// if(l==m)
						//	cout<<" "<<Jij_lattice_next[i][j]->pull_value()[l][m];
					}
					// cout<<endl;
					// cout<<endl;
				}
				Jij_lattice_next[i][j] = Jij_lattice_next[i][j]->pull_next();
				//	count++;
			}
			// cout<< "CONTEGGI "<< count <<endl;
			for (int l = 0; l < 3; l++)
				for (int m = 0; m < 3; m++)
					Jij_k[i][j][l][m] = C(Jij_tmp[l][m].real(), Jij_tmp[l][m].imag());
		}

	return Jij_k;
};

void lattice_J::free_fft(complex<double> ****Jij_k)
{
	int number_magnetic_atoms = spin_lattice->pull_number_magnetic_atoms();
	
	for (int i = 0; i < number_magnetic_atoms; i++)
		for (int j = 0; j < number_magnetic_atoms; j++)
			delete[] Jij_k[i][j];

	//for (int i = 0; i < number_atoms; i++)
		//delete[] Jij_k[i];
	
	delete[] Jij_k;
};

void lattice_J::print_fft(double *kpoint_tmp)
{
	int number_magnetic_atoms = spin_lattice->pull_number_magnetic_atoms();
	complex<double> ****Jij_k;
	Jij_k = pull_fft(kpoint_tmp);
	for (int i = 0; i < number_magnetic_atoms; i++)
		for (int j = 0; j < number_magnetic_atoms; j++)
		{
			cout << i << " " << j << endl;
			for (int l = 0; l < 3; l++)
			{
				for (int m = 0; m < 3; m++)
					cout << Jij_k[i][j][l][m].real() << "+i(" << Jij_k[i][j][l][m].imag() << ") ";
				cout << endl;
			}
		}

	free_fft(Jij_k);
};
///////// END DEFINITION CLASS J-LATTICE OF THE SYSTEM

///////// START DEFINITION CLASS HAMILTONIAN OF THE SYSTEM WITH RESPECT TO A PARTICULAR K POINT
class Hamiltonian_k
{
private:
	double *kpoint;
	/// the matrix shoud have 4 indices ijkl where i and j are associated to the creation+k and annihilation-k operator, and k and l are associated to the atoms of the crystal;
	/// for simplicity, we will consider a matrix kl with twice the dimension 2x2xkxl -> kxl
	complex<double> **H;
	lattice_spin *spin_crystal;
	lattice_J *J_crystal;
	int basis_dimension_H;
	//int flag_only_magneticatoms;
	/// Cholesky factorization Lower Triangular Matrix
	complex<double> **K;
	/// quantifying the gap opening to have H positive-deined (defined Cholesky factorization)
	double gap_opening;

public:
	Hamiltonian_k()
	{
		kpoint = nullptr;
		H = nullptr;
		spin_crystal = nullptr;
		J_crystal = nullptr;
		K = nullptr;
		basis_dimension_H = 0;
		//flag_only_magneticatoms = 0;
		gap_opening = 0;
	}
	void push_values(double *kpoint_tmp, lattice_J *J_crystal_tmp, lattice_spin *spin_crystal_tmp);
	void print();
	complex<double> **pull_Hk()
	{
		return H;
	}
	//void selecting_only_magneticatoms();
	// void Cholesky_decomposition();
	void Cholesky_decomposition_Lapack();
	// void SVD_decomposition();
	void checking_Cholesky();
	void Lapack_test();
	complex<double> **pull_Kk()
	{
		return K;
	}
	double *eigenvalues();
	complex<double> **eigenvectors();
	int pull_basis_dimension()
	{
		return basis_dimension_H;
	}
	void push_gap(double *epsilon)
	{
		for (int i = 0; i < basis_dimension_H; i++)
			H[i][i].real(H[i][i].real() + *epsilon);
		gap_opening = *epsilon;
	}
	//int pull_flag_only_magneticatoms()
	//{
	//	return flag_only_magneticatoms;
	//}
	void free_Hamiltonian_k();
};

void Hamiltonian_k::push_values(double *kpoint_tmp, lattice_J *J_crystal_tmp, lattice_spin *spin_crystal_tmp)
{

	spin_crystal = spin_crystal_tmp;
	J_crystal = J_crystal_tmp;
	int number_magnetic_atoms = spin_crystal->pull_number_magnetic_atoms();

	complex<double> factor2;
	factor2.real(2.0);
	factor2.imag(0.0);

	complex<double> A[number_magnetic_atoms][number_magnetic_atoms];
	complex<double> Aneg[number_magnetic_atoms][number_magnetic_atoms];
	complex<double> B[number_magnetic_atoms][number_magnetic_atoms];
	complex<double> C[number_magnetic_atoms][number_magnetic_atoms];

	int number_magnetic_atoms_2 = 2 * number_magnetic_atoms;
	basis_dimension_H = number_magnetic_atoms_2;
	// cout<<"fuori";

	if (H == nullptr)
	{
		// cout<<"dentro";
		H = new complex<double> *[number_magnetic_atoms_2];
		for (int i = 0; i < number_magnetic_atoms_2; i++)
			H[i] = new complex<double>[number_magnetic_atoms_2];
	}

	complex<double> ****Jij_kpos;
	complex<double> ****Jij_kneg;
	complex<double> ****Jij_knullptr;
	Jij_kpos = J_crystal->pull_fft(kpoint_tmp);

	double *kpoint_tmp2;
	kpoint_tmp2 = new double[3];
	double *kpoint_tmp3;
	kpoint_tmp3 = new double[3];
	for (int i = 0; i < 3; i++)
	{
		kpoint_tmp2[i] = (-1) * kpoint_tmp[i];
		kpoint_tmp3[i] = 0.0;
	}

	Jij_knullptr = J_crystal->pull_fft(kpoint_tmp3);
	Jij_kneg = J_crystal->pull_fft(kpoint_tmp2);
	delete[] (double *)kpoint_tmp2;
	delete[] (double *)kpoint_tmp3;

	for (int i=0;i<number_magnetic_atoms;i++)
	{
		for (int j=0;j<number_magnetic_atoms;j++)
		{
			complex<double> sum_1 = (0, 0);
			complex<double> sum_2 = (0, 0);
			complex<double> sum_5 = (0, 0);
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					sum_5 = sum_5 + spin_crystal->pull_spin()[i].pull_u()[k] * Jij_kpos[i][j][k][l] * conj(spin_crystal->pull_spin()[j].pull_u()[l]);
					sum_1 = sum_1 + spin_crystal->pull_spin()[i].pull_u()[k] * Jij_kneg[i][j][k][l] * conj(spin_crystal->pull_spin()[j].pull_u()[l]);
					sum_2 = sum_2 + spin_crystal->pull_spin()[i].pull_u()[k] * Jij_kneg[i][j][k][l] * spin_crystal->pull_spin()[j].pull_u()[l];
				}
			}
			double variable_tmp = sqrt(spin_crystal->pull_spin()[i].pull_modulus() * spin_crystal->pull_spin()[j].pull_modulus());
			Aneg[i][j] = (sum_5 / factor2) * variable_tmp;
			A[i][j] = (sum_1 / factor2) * variable_tmp;
			B[i][j] = (sum_2 / factor2) * variable_tmp;

			if (i != j)
			{
				C[i][j] = 0;
			}
			else
			{
				complex<double> sum_3 = (0, 0);
				complex<double> sum_4;
				for (int l = 0; l < number_magnetic_atoms; l++)
				{
					sum_4 = (0, 0);
					for (int k = 0; k < 3; k++)
						for (int r = 0; r < 3; r++)
							sum_4 = sum_4 + spin_crystal->pull_spin()[i].pull_orientation()[k] * Jij_knullptr[i][l][k][r] * spin_crystal->pull_spin()[l].pull_orientation()[r];
					sum_3 = sum_3 + sum_4 * spin_crystal->pull_spin()[l].pull_modulus();
				}
				C[i][j] = sum_3;
			}
		}
	}
	for (int i=0;i<number_magnetic_atoms_2;i++)
	{
		for (int j=0;j<number_magnetic_atoms_2;j++)
		{
			if ((i < number_magnetic_atoms) && (j < number_magnetic_atoms))
			{
				H[i][j] = factor2 * (A[i][j] - C[i][j]);
			}
			else if ((i < number_magnetic_atoms) && (j >= number_magnetic_atoms))
			{
				H[i][j] = factor2 * (B[i][j - number_magnetic_atoms]);
			}
			else if ((j < number_magnetic_atoms) && (i >= number_magnetic_atoms))
			{
				H[i][j] = factor2 * conj(B[j][i - number_magnetic_atoms]);
			}
			else
			{
				H[i][j] = factor2 * (Aneg[i - number_magnetic_atoms][j - number_magnetic_atoms] - C[i - number_magnetic_atoms][j - number_magnetic_atoms]);
			}
		}
	}

	J_crystal->free_fft(Jij_kneg);
	J_crystal->free_fft(Jij_kpos);
	J_crystal->free_fft(Jij_knullptr);
};

void Hamiltonian_k ::free_Hamiltonian_k()
{
	int n = basis_dimension_H;
	for (int i=0;i<n;i++)
	{
		delete[] (complex<double> *)H[i];
		delete[] (complex<double> *)K[i];
	}
	basis_dimension_H = 0;
	//flag_only_magneticatoms = 0;
	gap_opening = 0;
	H = nullptr;
	K = nullptr;
	kpoint = nullptr;
	spin_crystal = nullptr;
	J_crystal = nullptr;
};

void Hamiltonian_k::Cholesky_decomposition_Lapack()
{
	int n = basis_dimension_H;
	complex<double> **L_tot;
	L_tot = new complex<double> *[n];
	for (int i = 0; i < n; i++)
		L_tot[i] = new complex<double>[n];

	lapack_complex_double *H_tmp;
	H_tmp = (lapack_complex_double *)malloc(4 * n * n * sizeof(lapack_complex_double));

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			H_tmp[i * n + j] = H[i][j].real() + _Complex_I * (H[i][j].imag());

	// calculating the eigenvalues to check if the matrix is positive difined
	// this function is also calculating the Schur form
	char JOBVS = 'N';
	char SORT = 'N';
	LAPACK_Z_SELECT1 SELECT;
	int *SDIM;
	int N;
	int LDA;
	int LDVS;
	int matrix_layout = 101;
	int INFO;
	lapack_complex_double *W;
	lapack_complex_double *VS;
	VS = (lapack_complex_double *)malloc(n * n * sizeof(lapack_complex_double));
	W = (lapack_complex_double *)malloc(n * sizeof(lapack_complex_double));
	SDIM = (int *)malloc(sizeof(int));
	N = n;
	LDA = n;
	LDVS = n;

	INFO = LAPACKE_zgees(matrix_layout, JOBVS, SORT, SELECT, N, H_tmp, LDA, SDIM, W, VS, LDVS);	
	cout<<"Info pre-Cholesky "<<INFO<<endl;

	free(SDIM);
	free(VS);
	double positiveness = 1;
	for (int i = 0; i < n; i++)
		positiveness = positiveness * lapack_complex_double_real(W[i]);
	cout<<"determinant: "<<positiveness<<endl;
	int flag = 0;
	if (positiveness > 0)
	{
		cout<<"Positive defined "<<endl;
		flag = 1;
	}
	else if (positiveness == 0)
	{
		cout <<"Semi defined" << endl;
		flag = 2;
	}
	else
	{
		cout<<"Negative defined "<<endl;
		flag = 3;
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			H_tmp[i * n + j] = H[i][j].real() + _Complex_I * (H[i][j].imag());

	if (flag > 1)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				H_tmp[i * n + j] = (-1) * lapack_complex_double_real(H_tmp[i * n + j]) + _Complex_I * (lapack_complex_double_imag(H_tmp[i * n +j]));
	}

	char UPLO = 'U';
	INFO=LAPACKE_zpotrf2(matrix_layout, UPLO, N, H_tmp, LDA);
	cout<<"Info Cholesky "<<INFO<<endl;

	//double cycle_tmp=1;
	//int INFO2=INFO;
	//cout<<"Other attempt"<<endl;
	//while((INFO2>0)&&(cycle_tmp<=sqrt(INFO))){
	//	for (int i = 0; i < n; i++)
	//		for (int j = 0; j < n; j++){
	//			if((i<int(cycle_tmp))&&(j<int(cycle_tmp)))
	//				H_tmp[i * n + j] = (-1)*(H[i][j].real()) + _Complex_I * (H[i][j].imag());
	//			else
	//				H_tmp[i * n + j] = (H[i][j].real()) + _Complex_I * (H[i][j].imag());
	//		}
	//	INFO2=LAPACKE_zpotrf2(matrix_layout, UPLO, N, H_tmp, LDA);
	//	cout<<"Info Cholesky "<<INFO2<<endl;
	//	cycle_tmp=cycle_tmp+1;
	//}
	

	typedef std::complex<double> C;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (j >= i)
			{
				L_tot[i][j].real(lapack_complex_double_real(H_tmp[i * n + j]));
				L_tot[i][j].imag(lapack_complex_double_imag(H_tmp[i * n + j]));
			}
			else
				L_tot[i][j] = C(0, 0);
		}
	complex<double> L_tot_T[n][n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			L_tot_T[i][j] = L_tot[j][i];
	
	free(H_tmp);
	free(W);
	K = L_tot;
};

void Hamiltonian_k::Lapack_test()
{
	/// Simple conjugation to check the working of lapack library
	int n = basis_dimension_H;
	lapack_complex_float *vector = new lapack_complex_float[n];
	lapack_int *N = new lapack_int[1];
	lapack_int flag = 1;
	lapack_int IDIST;
	lapack_int ISEED;
	IDIST = 1;
	*N = n;
	ISEED = 2;

	cout << "FLAG: " << flag << endl;
	for (int j = 0; j < n; j++)
	{
		cout << lapack_complex_float_real(vector[j]) << "+i(" << lapack_complex_float_imag(vector[j]) << ") ";
	}
	cout << endl;
	flag = LAPACKE_clacgv(*N, vector, 1);
	for (int j = 0; j < n; j++)
		cout << lapack_complex_float_real(vector[j]) << "+i(" << lapack_complex_float_imag(vector[j]) << ") ";
	cout << endl;

	cout << "FLAG: " << flag << endl;

	delete[] vector;
};

void Hamiltonian_k::checking_Cholesky()
{
	int n = basis_dimension_H;
	complex<double> M[n][n];
	typedef std::complex<double> C;
	cout << "Testing Cholesky decomposition " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			M[i][j] = C(0, 0);
			for (int k = 0; k < n; k++)
				M[i][j] = M[i][j] + conj(K[k][i]) * (K[k][j]);
			cout << M[i][j].real() << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << "K matrix from the Cholesky decomposition " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << K[i][j] << " ";
		}
		cout << endl;
	}
};

double *Hamiltonian_k::eigenvalues()
{
	try
	{
		if (K == nullptr)
			throw 1;
	}
	catch (int tmp)
	{
		throw std::system_error{errno, std::generic_category(), "failed Cholesky"};
	}

	int n = basis_dimension_H;
	int n_2 = basis_dimension_H / 2;

	complex<double> **bosonic_H;
	complex<double> **g;
	bosonic_H = new complex<double> *[n];
	g = new complex<double> *[n];
	for (int i = 0; i < n; i++)
	{
		bosonic_H[i] = new complex<double>[n];
		g[i] = new complex<double>[n];
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			g[i][j].imag(0.0);
			bosonic_H[i][j].real(0.0);
			bosonic_H[i][j].imag(0.0);
			if (i == j)
			{
				if ((i < n_2) && (j < n_2))
					g[i][j].real(1.0);
				else
					g[i][j].real(-1.0);
			}
			else
				g[i][j].real(0.0);
		}

	/// applying the Colpa method to preserve bosonic commutation operations
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int l = 0; l < n; l++)
				for (int m = 0; m < n; m++)
					bosonic_H[i][j] = bosonic_H[i][j] + K[i][l] * g[l][m] * conj(K[j][m]);

	//for (int i = 0;i<n; i++){
	//	for (int j = 0;j<n; j++)
	//		cout<<bosonic_H[i][j]<<" ";
	//	cout<<endl;
	//}

	lapack_complex_double *temporary;
	temporary = (lapack_complex_double*) malloc(n*n*sizeof(lapack_complex_double));

	for (int i = 0;i<n; i++){
		for (int j = 0;j<n; j++)
			temporary[i * n + j]= bosonic_H[i][j].real()+ _Complex_I*(bosonic_H[i][j].imag());
	}

	for (int i = 0; i < n; i++)
		delete[] (complex<double> *)bosonic_H[i];
	delete[] (complex<double> **)bosonic_H;
	

	int N=n;
	int LDA=n;
	int matrix_layout=101;
	int INFO;
	lapack_complex_double *W; 
	char JOBVL = 'N'; char JOBVR = 'N';
	lapack_complex_double *VL; lapack_complex_double *VR;
	int LDVL=1; int LDVR=1;
	W = (lapack_complex_double *)malloc(n * sizeof(lapack_complex_double));

	INFO=LAPACKE_zgeev(matrix_layout,JOBVL,JOBVR,N,temporary,LDA,W,VL,LDVL,VR,LDVR);
	
	cout<<"RESULTS"<<endl;
	for (int i = 0; i < n; i++)
		cout<<lapack_complex_double_real(W[i])<<" ";
	cout<<endl;

	int check=0;
	for(int i=0;i<n;i++){
		//cout<<lapack_complex_double_real(W[i])<<" "<<lapack_complex_double_imag(W[i])<<" ";
		if(lapack_complex_double_real(W[i])<0){
			W[i]=(-1)*lapack_complex_double_real(W[i])+ _Complex_I*(lapack_complex_double_imag(W[i]));
			check++;
		}	
	}

	cout<<endl;

	if (check>n_2){
		cout<<"negative definite hamiltonian"<<endl;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				temporary[i * n + j]=(-1)*lapack_complex_double_real(temporary[i * n + j])+_Complex_I*(lapack_complex_double_imag(temporary[i * n + j]));
		INFO=LAPACKE_zgeev(matrix_layout,JOBVL,JOBVR,N,temporary,LDA,W,VL,LDVL,VR,LDVR);
		check=0;
		for(int i=0;i<n;i++)
			if(lapack_complex_double_real(W[i])<0){
				W[i]=(-1)*lapack_complex_double_real(W[i])+ _Complex_I*(lapack_complex_double_imag(W[i]));
				check++;
			}
		if (check>n_2)
			cout<<"problem diagonalization"<<endl;
	}else if(check<n_2)
		cout<<"increase epsilon"<<endl;
	else
		cout<<"positive definite hamiltonian"<<endl;

	delete[] (lapack_complex_double *)temporary;
	double* WR=new double[n];	
	for(int i=0;i<n;i++)
		WR[i]=lapack_complex_double_real(W[i]);	
	
	mergeSort(WR, 0, n - 1);
	
	double *w_final;
	w_final = new double[n_2];
	int counting=0;
	for (int i = 0; i < n; i++){
		if(i%2==0){
			w_final[counting] = WR[i];
			counting++;
		}
	}
	free(W);
	free(WR);
	for (int i = 0; i < n; i++)
		delete[] (complex<double> *)g[i];
	delete[] (complex<double> **)g;
	return w_final;
};

complex<double> **Hamiltonian_k::eigenvectors()
{
	try
	{
		if (K == nullptr)
			throw 1;
	}
	catch (int tmp)
	{
		throw std::system_error{errno, std::generic_category(), "failed Cholesky"};
	}

	int n = basis_dimension_H;
	int n_2 = basis_dimension_H / 2;

	complex<double> **bosonic_H;
	complex<double> **g;
	bosonic_H = new complex<double> *[n];
	g = new complex<double> *[n];
	for (int i = 0; i < n; i++)
	{
		bosonic_H[i] = new complex<double>[n];
		g[i] = new complex<double>[n];
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			g[i][j].imag(0.0);
			bosonic_H[i][j].real(0.0);
			bosonic_H[i][j].imag(0.0);
			if (i == j)
			{
				if ((i < n_2) && (j < n_2))
					g[i][j].real(1.0);
				else
					g[i][j].real(-1.0);
			}
			else
				g[i][j].real(0.0);
		}

	/// applying the Colpa method to preserve bosonic commutation operations
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int l = 0; l < n; l++)
				for (int m = 0; m < n; m++)
					bosonic_H[i][j] = bosonic_H[i][j] + K[i][l] * g[l][m] * conj(K[j][m]);

	lapack_complex_double *temporary;
	temporary = (lapack_complex_double *)malloc(n * n * sizeof(lapack_complex_double));

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			temporary[i * n + j] = bosonic_H[i][j].real() + _Complex_I * (H[i][j].imag());

	for (int i = 0; i < n; i++)
		delete[] (complex<double> *)bosonic_H[i];
	delete[] (complex<double> **)bosonic_H;
	for (int i = 0; i < n; i++)
		delete[] (complex<double> *)g[i];
	delete[] (complex<double> **)g;

	int N;
	int LDA;
	int matrix_layout = 101;
	int INFO;
	/// WR are the eigenvalues
	// U are the eigenvectors
	double *WR;
	complex<double> **U;
	char JOBZ = 'V';
	char UPLO = 'L';
	// saving all the eigenvalues
	WR = (double *)malloc(n * sizeof(double));

	// I am considering the fact that the combination could include a phase
	N = n;
	LDA = n;
	INFO = LAPACKE_zheev(matrix_layout, JOBZ, UPLO, N, temporary, LDA, WR);

	complex<double> **U_final;
	U_final = new complex<double> *[n_2];
	for (int i = 0; i < n; i++)
		U_final[i] = new complex<double>[n_2];

	double WR_tmp[n][2];
	for (int i = 0; i < n; i++)
	{
		WR_tmp[i][0] = WR[i];
		WR_tmp[i][1] = 0;
		// cout<<WR[i]<<endl;
	}

	/// ordering eigenvalues to make ordering of the associated eigenvector consistent during the path for the Berry curvature
	mergeSort(WR, 0, n - 1);

	/// considering only the first (negative) eigenvalues
	int ordering[n_2];

	// cout<<"START"<<endl;
	int flag;
	int j;
	for (int i = 0; i < n_2; i++)
	{
		j = 0;
		flag = 0;
		while (flag != 1)
		{
			// cout<<WR[i]<<" "<<WR_tmp[j][0]<<endl;
			if ((WR_tmp[j][0] == WR[i]) && (WR_tmp[j][1] != 1))
			{
				ordering[i] = j;
				WR_tmp[j][1] = 1;
				flag = 1;
				// cout<<i<<" "<<j<<" "<<ordering[i]<<endl;
			}
			else
				j++;
		}
	}
	// cout<<"END"<<endl;

	for (int i = 0; i < n_2; i++)
	{
		for (int j = 0; j < n_2; j++)
		{
			U_final[i][j] = temporary[j * n + ordering[i]];
			// cout<<U_final[i][j]<<" ";
		}
		// cout<<endl;
	}
	double norm_tmp;
	complex<double> norm;
	norm.imag(0.0);
	for(int i = 0; i< n_2; i++){
		norm=0;
		for(int j = 0; j< n_2; j++){
			norm_tmp=norm.real()+abs(U_final[i][j]);
			norm.real(norm_tmp);
		}
		norm.real(sqrt(norm.real()));
		for(int j = 0; j< n_2; j++)
			if(norm.real()!=0)
				U_final[i][j]=U_final[i][j]/norm;
	}

	free(WR);
	free(temporary);

	return U_final;
};


void Hamiltonian_k::print()
{
	cout << endl;
	cout << "Real Part Hamiltonian" << endl;
	for (int i = 0; i < basis_dimension_H; i++)
	{
		for (int j = 0; j < basis_dimension_H; j++)
			printf("%.4f ", H[i][j].real());
		cout << endl;
	}
	cout << endl;
	cout << "Immaginary Part Hamiltonian" << endl;
	for (int i = 0; i < basis_dimension_H; i++)
	{
		for (int j = 0; j < basis_dimension_H; j++)
			printf("%.4f ", H[i][j].imag());
		cout << endl;
	}
};
///////// END DEFINITION CLASS HAMILTONIAN OF THE SYSTEM WITH RESPECT TO A PARTICULAR K POINT


///////// START DEFINITION CLASS K_POINTS: LIST OF K POINTS ALONG THE PATH POINTED OUT BY THE FILE K_POINTS 
class K_points
{
private:
	int number_special_k_points;
	double **special_k_points;
	int *steps_k_points;
	double **list_k_points;
	int total_number_k_points;

public:
	K_points()
	{
		number_special_k_points = 0;
		special_k_points = nullptr;
		list_k_points = nullptr;
		steps_k_points = nullptr;
		total_number_k_points = 0;
	}
	void push_values(ifstream *k_points_file_tmp);
	int pull_total_number_k_points()
	{
		return total_number_k_points;
	}
	double **pull_list_k_points()
	{
		return list_k_points;
	}
};

void K_points::push_values(ifstream *k_points_file_tmp)
{
	try
	{
		if (k_points_file_tmp == nullptr)
			throw 1;
	}
	catch (int flag)
	{
		throw std::system_error{errno, std::generic_category(), "k_points_file_error"};
	}
	/// here I read the file twice; in order to avoid improper complexity into the code
	/// being this file relativily short
	string line_tmp;
	k_points_file_tmp->seekg(0);
	while (k_points_file_tmp->peek() != EOF)
	{
		getline(*k_points_file_tmp, line_tmp);
		number_special_k_points++;
	}
	special_k_points = new double *[number_special_k_points];
	steps_k_points = new int[number_special_k_points];
	for (int i = 0; i < number_special_k_points; i++)
		special_k_points[i] = new double[3];

	int count = 0;
	k_points_file_tmp->seekg(0);
	for (int m = 0; m < number_special_k_points; m++)
	{
		for (int l = 0; l < 3; l++)
			*k_points_file_tmp >> special_k_points[m][l];
		*k_points_file_tmp >> steps_k_points[m];
		total_number_k_points = total_number_k_points + steps_k_points[m];
	}

	total_number_k_points++;
	list_k_points = new double *[total_number_k_points];
	for (int l = 0; l < total_number_k_points; l++)
		list_k_points[l] = new double[3];

	count = 0;
	for (int i = 0; i < number_special_k_points; i++)
	{
		for (int j = 0; j < steps_k_points[i]; j++)
			for (int l = 0; l < 3; l++)
				list_k_points[count + j][l] = (special_k_points[i + 1][l] - special_k_points[i][l]) * ((double)j / steps_k_points[i]) + special_k_points[i][l];
		count = count + steps_k_points[i];
	}
	list_k_points[count] = special_k_points[number_special_k_points - 1];
};
///////// END DEFINITION CLASS K_POINTS: LIST OF K POINTS ALONG THE PATH POINTED OUT BY THE FILE K_POINTS 

///////// START DEFINITION CLASS K_POINT GRID USED TO CALCULATE BERRY CURVATURE
class K_points_grid
{
private:
	int plane[3];
	double spacing;
	int number_k_points[2];
	double *shift;
	double ***k_point_grid;
	double **vector_reciprocal;

public:
	K_points_grid()
	{
		spacing = 0;
		vector_reciprocal = nullptr;
		k_point_grid = nullptr;
		shift = new double[3];
		for (int i = 0; i < 3; i++)
		{
			shift[i] = 0.0;
			plane[i] = 0;
		}
	}
	void K_points_grid_push_values(lattice_crystal *lattice_crystal_tmp, double spacing_tmp, double *shift_tmp, int *plane_tmp);
	double *K_points_grid_get_values_grid(int i_tmp, int j_tmp)
	{
		return k_point_grid[i_tmp][j_tmp];
	}
	double K_points_grid_get_vaues_spacing()
	{
		return spacing;
	}
	double *K_points_grid_get_vaues_shift()
	{
		return shift;
	}
	int *K_points_grid_get_number_k_points()
	{
		return number_k_points;
	}
	int *K_points_grid_get_plane()
	{
		return plane;
	}
	double **K_points_grid_get_vector_reciprocal()
	{
		return vector_reciprocal;
	}
	void K_points_grid_print();
};

void K_points_grid::K_points_grid_push_values(lattice_crystal *lattice_crystal_tmp, double spacing_tmp, double *shift_tmp, int *plane_tmp)
{
	vector_reciprocal = new double *[2];
	double factor = 2 * pigreco / lattice_crystal_tmp->pull_volume();
	spacing = spacing_tmp;
	shift = shift_tmp;

	int count = 0;
	for (int i = 0; i < 3; i++)
	{
		if ((plane_tmp[i] != 0) && (count < 2))
		{
			if (i == 0)
			{
				// cout<<"1"<<endl;
				vector_reciprocal[count] = function_vector_product(lattice_crystal_tmp->pull_bravais_lattice()[1], lattice_crystal_tmp->pull_bravais_lattice()[2]);
			}
			else if (i == 1)
			{
				// cout<<"2"<<endl;
				vector_reciprocal[count] = function_vector_product(lattice_crystal_tmp->pull_bravais_lattice()[2], lattice_crystal_tmp->pull_bravais_lattice()[0]);
			}
			else
			{
				// cout<<"3"<<endl;
				vector_reciprocal[count] = function_vector_product(lattice_crystal_tmp->pull_bravais_lattice()[0], lattice_crystal_tmp->pull_bravais_lattice()[1]);
			}
			count++;
		}
	}

	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 3; j++)
			vector_reciprocal[i][j] = factor * vector_reciprocal[i][j];

	/// while the kpoints are in units of 1/a.u. the shift is in units of the spacing
	number_k_points[0] = (function_module(vector_reciprocal[0]) / spacing);
	number_k_points[1] = (function_module(vector_reciprocal[1]) / spacing);

	// cout<<number_k_points[0]<<" "<<number_k_points[1]<<endl;

	k_point_grid = new double **[number_k_points[0]];
	for (int i = 0; i < number_k_points[0]; i++)
	{
		k_point_grid[i] = new double *[number_k_points[1]];
		for (int j = 0; j < number_k_points[1]; j++)
			k_point_grid[i][j] = new double[3];
	}

	for (int i = 0; i < number_k_points[0]; i++)
		for (int j = 0; j < number_k_points[1]; j++)
			for (int r = 0; r < 3; r++)
				k_point_grid[i][j][r] = ((double)i / number_k_points[0]) * (shift[r] + vector_reciprocal[0][r]) + ((double)j / number_k_points[1]) * (shift[r] + vector_reciprocal[1][r]);
};

void K_points_grid::K_points_grid_print()
{
	for (int i = 0; i < number_k_points[0]; i++)
	{
		for (int j = 0; j < number_k_points[1]; j++)
		{
			cout << " ( ";
			for (int r = 0; r < 3; r++)
				cout << k_point_grid[i][j][r] << " ";
			// cout<<i<<" "<<j;
			cout << " ) ";
		}
		cout << endl;
	}
};
///////// END DEFINITION CLASS K_POINT GRID USED TO CALCULATE BERRY CURVATURE

///////// START DEFINITION CLASS BERRY CURVATURE
class Berry_curvature
{
private:
	K_points_grid *K_grid;
	int O_basis_dimension;
	double *epsilon;
	lattice_spin *Crystal_spin;
	lattice_J *Crystal_J;
	complex<double> ***Omega;

public:
	Berry_curvature()
	{
		K_grid = nullptr;
		O_basis_dimension = 0;
		Crystal_spin = nullptr;
		Crystal_J = nullptr;
		Omega = nullptr;
	};
	void Berry_curvature_push_values(K_points_grid *K_grid_tmp, lattice_J *Crystal_J_tmp, lattice_spin *Crystal_spin_tmp, double *epsilon_tmp, int band_number);
	void Berry_curvature_print(ofstream *Re_Omega_file, ofstream *Im_Omega_file);
	complex<double> Berry_curvature_k(double *k_point,int number_band);
};

void Berry_curvature ::Berry_curvature_push_values(K_points_grid *K_grid_tmp, lattice_J *Crystal_J_tmp, lattice_spin *Crystal_spin_tmp, double *epsilon_tmp, int band_number)
{
	K_grid = K_grid_tmp;
	Crystal_J = Crystal_J_tmp;
	Crystal_spin = Crystal_spin_tmp;
	epsilon = epsilon_tmp;
	/// considering only magnetic atoms
	O_basis_dimension = Crystal_spin->pull_number_magnetic_atoms() * 2;
	int *number_k_points = K_grid->K_points_grid_get_number_k_points();

	int n_2 = O_basis_dimension / 2;
	int n = O_basis_dimension;
	int k0 = number_k_points[0];
	int k1 = number_k_points[1];

	double *k_point = new double[3];

	complex<double> ****uk;
	uk = new complex<double>*** [k0];
	for (int i = 0; i < k0; i++){
		uk[i]=new complex<double>**[k1];
	}
	
	Hamiltonian_k *Hk;
	Hk = new Hamiltonian_k[1];
	
	for (int i = 0; i < k0; i++)
		for (int j = 0; j < k1; j++){
			Hk->push_values(K_grid->K_points_grid_get_values_grid(i,j), Crystal_J, Crystal_spin);
			Hk->push_gap(epsilon);
			Hk->Cholesky_decomposition_Lapack();
			uk[i][j] = Hk->eigenvectors();
			Hk->free_Hamiltonian_k();
		}

	Omega = new complex<double> **[k0];
	for (int i = 0; i < k0; i++){
		Omega[i] = new complex<double>*[k1];
		for (int j = 0; j < k1; j++)
			Omega[i][j]=new complex<double>[n_2];
	}

	complex<double> ***uk_tmp;
	uk_tmp = new complex<double>**[4];

	complex<double> **phase1;
	phase1 = new complex<double> *[n_2];
	complex<double> **phase2;
	phase2 = new complex<double> *[n_2];
	complex<double> **phase3;
	phase3 = new complex<double> *[n_2];
	complex<double> **phase4;
	phase4 = new complex<double> *[n_2];
	for (int i = 0; i < n_2; i++)
	{
		phase1[i] = new complex<double>[n_2];
		phase2[i] = new complex<double>[n_2];
		phase3[i] = new complex<double>[n_2];
		phase4[i] = new complex<double>[n_2];
	}

	complex<double> phase1_tmp;
	complex<double> phase2_tmp;
	complex<double> phase3_tmp;
	complex<double> phase4_tmp;

	for (int i = 0; i < k0; i++)
		for (int j = 0; j < k1; j++)
		{
			cout << i << " / " << k0 << endl;
			// cout<<j<<" / "<<k1<<endl;
			//		////periodic boundary conditions
			if ((i < number_k_points[1] - 1) && (j < number_k_points[0] - 1))
			{
				uk_tmp[0] = uk[i][j];
				uk_tmp[3] = uk[i][j + 1];
				uk_tmp[2] = uk[i + 1][j+1];
				uk_tmp[1] = uk[i + 1][j];
			}
			else if ((i < number_k_points[1] - 1) && (j == number_k_points[1] - 1))
			{
				uk_tmp[0] = uk[i][j];
				uk_tmp[3] = uk[i][0];
				uk_tmp[2] = uk[i + 1][0];
				uk_tmp[1] = uk[i + 1][j];
			}
			else if ((i == number_k_points[1] - 1) && (j < number_k_points[1] - 1))
			{
				uk_tmp[0] = uk[i][j];
				uk_tmp[3] = uk[i][j + 1];
				uk_tmp[2] = uk[0][j + 1];
				uk_tmp[1] = uk[0][j];
			}
			else
			{
				uk_tmp[0] = uk[i][j];
				uk_tmp[3] = uk[i][0];
				uk_tmp[2] = uk[0][0];
				uk_tmp[1] = uk[0][j];
			}
			for (int l = 0; l < n_2; l++)
				for (int m = 0; m < n_2; m++)
				{
					phase1[l][m] = function_scalar_product_2(uk_tmp[0], uk_tmp[1], l, m, n_2);
					phase2[l][m] = function_scalar_product_2(uk_tmp[1], uk_tmp[2], l, m, n_2);
					phase3[l][m] = function_scalar_product_2(uk_tmp[3], uk_tmp[2], l, m, n_2);
					phase4[l][m] = function_scalar_product_2(uk_tmp[3], uk_tmp[0], l, m, n_2);
				}
			
			for(int l=0;l<n_2;l++){
				phase1_tmp = function_quantum_of_phase(arg(phase1[i][i]));
				phase2_tmp = function_quantum_of_phase(arg(phase2[i][i]));
				phase3_tmp = function_quantum_of_phase(arg(phase3[i][i]));
				phase4_tmp = function_quantum_of_phase(arg(phase4[i][i]));
				Omega[i][j][l] = phase1_tmp + phase2_tmp - phase3_tmp - phase4_tmp;
			}
		}

	for(int r=0;r<4;r++)
		delete[](complex<double>*)uk[r];
	delete[] (complex<double> **)uk;

	for (int r = 0; r < n_2; r++)
	{
		delete[] (complex<double> *)phase1[r];
		delete[] (complex<double> *)phase2[r];
		delete[] (complex<double> *)phase3[r];
		delete[] (complex<double> *)phase4[r];
	}
	delete[] (complex<double> **)phase1;
	delete[] (complex<double> **)phase2;
	delete[] (complex<double> **)phase3;
	delete[] (complex<double> **)phase4;
	delete[] Hk;
};

complex<double> Berry_curvature ::Berry_curvature_k(double *k_point,int number_band)
{
	////finding the k points nearer to the input k_point and doing an interpolation
	///(the interpolation is done considering only a part of the k points around, because of the laziness of the programmer to consider all the cases)
	/// the main tep is to individuate the point in the Kpoint_grid
	int count = 0;
	double k_point_tmp[2];
	while (count < 2)
	{
		for (int i = 0; i < 3; i++)
			k_point_tmp[count] = k_point[i] * K_grid->K_points_grid_get_vector_reciprocal()[count][i];
		count++;
	}
	int i = (int)(k_point_tmp[0] / K_grid->K_points_grid_get_vaues_spacing());
	int j = (int)(k_point_tmp[1] / K_grid->K_points_grid_get_vaues_spacing());
	int k0 = K_grid->K_points_grid_get_number_k_points()[0];
	int k1 = K_grid->K_points_grid_get_number_k_points()[1];

	return Omega[i][j][number_band];
};

void Berry_curvature ::Berry_curvature_print(ofstream *Re_Omega_file, ofstream *Im_Omega_file)
{
	int *number_k_points = K_grid->K_points_grid_get_number_k_points();
	for (int i = 0; i < number_k_points[0]; i++)
	{
		for (int j = 0; j < number_k_points[1]; j++)
		{
			for (int r=0; r< O_basis_dimension/2;r++){
				(*Re_Omega_file) << i << " " << j << " " << Omega[i][j][r].real()<<endl;
				(*Im_Omega_file) << i << " " << j << " " << Omega[i][j][r].imag()<<endl;
			}
		}
		(*Re_Omega_file) << endl;
		(*Im_Omega_file) << endl;
	}
};
///////// END DEFINITION CLASS BERRY CURVATURE

int main(){
	////START CONSTRUCTION CRYSTAL
	////atoms coordinates are in the same units of the bravais lattice vectors
	/////// remind to write also the files spin.txt coupling.txt atoms.txt and so on... before running the program
	///// inizializzazione di questo oggetto richiede un file structure.txt ed un file basis.txt
	///// structure.txt = reticolo di Bravais | basis.txt == coordinate atomi (ordine consistente con la nomenclatura atomi del file coupling.txt)
	lattice_crystal crystal;
	ifstream bravais_lattice_file;
	bravais_lattice_file.open("bravais.lattice.data");
	ifstream atoms_coordinates_file;
	atoms_coordinates_file.open("atoms.coordinates.data");
	crystal.push_values(&bravais_lattice_file, &atoms_coordinates_file);
	bravais_lattice_file.close();
	atoms_coordinates_file.close();
	//crystal.print();
	///// inizializzazione di questo oggetto richiede un file spin.txt
	///// spin.txt == spin atomi (ordine consistente con la nomenclatura atomi del file coupling.txt)
	lattice_spin crystal_spin;
	ifstream atoms_spins_file;
	atoms_spins_file.open("spins.data");
	crystal_spin.push_values(&crystal, &atoms_spins_file);
	atoms_spins_file.close();
	//crystal_spin.print();
	int number_atoms=crystal.pull_number_atoms();
	for(int i=0;i<number_atoms;i++)
		cout<<i<<" "<<crystal_spin.pull_ordering_from_crystal_to_magnetic()[i]<<endl;
	int number_magnetic_atoms=crystal_spin.pull_number_magnetic_atoms();
	for(int i=0;i<number_magnetic_atoms;i++)
		cout<<i<<" "<<crystal_spin.pull_ordering_from_magnetic_to_crystal()[i]<<endl;
	//int number_magnetic_atoms=crystal_spin.pull_number_magnetic_atoms();
	///// inizializzazione di questo oggetto richiede un file coupling.txt
	///// coupling.txt ==  elemento 1 (elementi :qindicati con un indice 0,1,2,...)| elemento 2 | distanza coupling (in angstrom) | J_value(isotropic)
	double max_distance_coupling = 5.500;
	lattice_J crystal_J;
	ifstream couplings_file;
	couplings_file.open("couplings.data");
	crystal_J.push_values(max_distance_coupling,&crystal_spin,&couplings_file);
	couplings_file.close();
	//crystal_J.print();
	//crystal_J.building_values(max_distance_coupling,&crystal_spin);
	//crystal_J.print();
	////END CONSTRUCTION CRYSTAL

	//////START STUDYING THE DIAGONALIZATION FOR ONE K POINT (000)
	////////THIS PART IS NEEDED TO DEFINE THE BASIS OF ONLY MAGNETIC ATOMS (FOR BERRY PHASE)
	double *k_point;
	k_point = new double[3];
	k_point[0] = 0.0;
	k_point[1] = 1.0;
	k_point[2] = 0.0;
	////////crystal_J.print_fft(k_point);
	Hamiltonian_k Hk;
	Hk.push_values(k_point, &crystal_J, &crystal_spin);
	Hk.print();
	//Hk.selecting_only_magneticatoms();
	double *epsilon = new double;
	///2.7 needed in order to avoid immaginary solutions
	*epsilon = 0.1;
	Hk.push_gap(epsilon);
	Hk.print();
	Hk.Cholesky_decomposition_Lapack();
	Hk.checking_Cholesky();
	cout<<"Magnons for k: "<<k_point[0]<< k_point[1]<< k_point[2]<<endl;
	double* w_k;
	w_k=Hk.eigenvalues();
	for(int i=0;i<number_magnetic_atoms;i++)
		cout<<w_k[i]<<" ";
	cout<<endl;
	////complex<double> **uk;
	////uk = Hk.eigenvectors();
	////for (int j = 0; j < number_magnetic_atoms; j++)
	////{
	////	for (int i = 0; i < number_magnetic_atoms; i++)
	////		cout << uk[i][j] << " ";
	////	cout << endl;
	////}
	////////END STUDYING THE DIAGONALIZATION FOR ONE K POINT

	//////TEST DETERMINANT
	// cout<<"TEST DETERMINANT"<<endl;
	// int r=Hk.pull_basis_dimension();
	// for(int i=0;i<r;i++){
	//	for(int j=0;j<r;j++)
	//		cout<<Hk.pull_Hk()[i][j]<<" ";
	//	cout<<endl;
	// }
	// cout<<endl;
	// cout<<function_determinant(Hk.pull_Hk(),r);

	////START STUDYING THE DIAGONALIZATION FOR A PATH OF K POINTS
	/// file esterno kpoints.txt : special K points | number of points between
	// K1 4
	// K2 5
	/// K3 0 	k1,4points..,K2,5points,K3
	/// initializing k points list of the path
	ifstream k_points_file_tmp;
	k_points_file_tmp.open("kpoints.data");
	K_points k_points;
	k_points.push_values(&k_points_file_tmp);
	k_points_file_tmp.close();
	double **list_k_points;
	list_k_points = k_points.pull_list_k_points();
	int number_k_points;
	number_k_points = k_points.pull_total_number_k_points();
	//cout<<"List of K points: "<<endl;
	//for(int i=0;i<number_k_points;i++){
	//	for(int r=0;r<3;r++)
	//		cout<<list_k_points[i][r]<<" ";
	//	cout<<endl;
	//}
	//// calculating bands along the k path underlined in the file k_points_file
	// saving the bands in the file magnons.data
	double **w;
	w = new double *[number_k_points];
	ofstream magnons_file_tmp;
	magnons_file_tmp.open("magnons.data");
	Hamiltonian_k *H;
	H = new Hamiltonian_k[number_k_points];
	complex<double> ***u;
	u = new complex<double> **[number_k_points];
	//// parallelizing on the k calculations
	//#pragma omp parallel for
	for (int i = 0; i < number_k_points; i++)
	{
		w[i] = new double[number_magnetic_atoms];
		H[i].push_values(list_k_points[i], &crystal_J, &crystal_spin);
		H[i].push_gap(epsilon);
		H[i].Cholesky_decomposition_Lapack();
		w[i] = H[i].eigenvalues();
	}
	//saving the bands in the file magnons.data
	for (int i = 0; i < number_k_points; i++)
	{
		for (int m = 0; m < number_magnetic_atoms; m++)
			magnons_file_tmp << w[i][m] << "	";
		magnons_file_tmp << endl;
	}
	magnons_file_tmp.close();
	// cout<<"Eigenvectors: "<<endl;
	// for(int k=0;k<number_k_points;k++){
		//cout<<"K point: "<<k<<endl;
		//for(int j=0;j<number_magnetic_atoms*2;j++){
		//	for(int i=0;i<number_magnetic_atoms;i++)
		//		cout<<u[k][i][j]<<" ";
		//cout<<endl;
		//}
		//cout<<endl;
	// }
	//END STUDYING THE DIAGONALIZATION FOR A PATH OF K POINTS

	//////START STUDYING THE BERRY CURVATURE IN A K PLANE
	//K_points_grid K_grid;
	//double spacing_grid = 0.01;
	//double *shift_grid;
	//shift_grid = new double[3];
	//for (int i = 0; i < 3; i++)
	//	shift_grid[i] = 0.0;
//
	//int *plane_grid;
	//int band_number = 1;
	//plane_grid = new int[3];
	///// plane xy
	//plane_grid[0] = 1;
	//plane_grid[1] = 1;
	//plane_grid[2] = 0;
	/// creating the K grid in xy plane
	//K_grid.K_points_grid_push_values(&crystal, spacing_grid, shift_grid, plane_grid);
	//K_grid.K_points_grid_print();
	//Berry_curvature Omega;
	//Omega.Berry_curvature_push_values(&K_grid, &crystal_J, &crystal_spin, epsilon, band_number);
	//ofstream Re_Omega_file;
	//ofstream Im_Omega_file;
	//Re_Omega_file.open("Re_Berry_curvature.data");
	//Im_Omega_file.open("Im_Berry_curvature.data");
	//Omega.Berry_curvature_print(&Re_Omega_file, &Im_Omega_file);
	//Im_Omega_file.close();
	//Re_Omega_file.close();
	//// considering the same k path of the bands...
	//// complex<double>* Tr_Omega_k;
	//// Tr_Omega_k=new complex<double>;
	//// Tr_Omega_k=Omega.Berry_curvature_k(k_point);
	///// END STUDYING THE BERRY CURVATURE IN A K PLANE
//
	///// START STUDYING THE BERRY CURVATURE ALONG THE K PATH OF THE MAGNONS SPECTRA
	///// file Berry_path.data
	//double **O;
	//O = new double *[number_k_points];
	//ofstream Berry_file_tmp;
	//Berry_file_tmp.open("Berry_path.data");
	//complex<double> value_tmp;
//
	//for (int i = 0; i < number_k_points; i++)
	//{
	//	value_tmp = Omega.Berry_curvature_k(list_k_points[i]);
	//	Berry_file_tmp << value_tmp.imag() << endl;
	//}
//
	//Berry_file_tmp.close();
	///// END STUDYING THE BERRY CURVATURE ALONG THE K PATH OF THE MAGNONS SPECTRA
	return 0;
}

