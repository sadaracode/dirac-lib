#pragma once

#ifndef __UTILS_H
#define __UTILS_H

#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

#include <mkl.h>
#include <complex>
#include "dirac.h"
#include <intrin.h>
#include <stdlib.h>
#define Ne 18
#define LDA Ne
#define LDZ Ne

#define PI 3.14159265

using namespace std;

template <typename T> class Operator;
template <typename T> class State;

//typedef double realnum;

/* Complex datatype */
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;

#pragma endregion

namespace Utils {
	enum class Qtype {
		D,
		X,
		Z,
		O
	};

	template<typename T>
	Operator<T> annihilation(const T& n);

	template<typename T>
	Operator<T> creation(const T& n);

	template<typename T>
	Operator<T> pauliI(const T& n);

	template<typename T>
	Operator<T> pauliI();

	template<typename T>
	Operator<T> pauliX();

	template<typename T>
	Operator<T> pauliZ();

	template<typename T>
	Operator<T> pauliY();

	template<typename T>
	Operator<T> pauliPlus();

	template<typename T>
	Operator<T> pauliMinus();

	template<typename T>
	Operator<T> Hadamard();

	template<typename T>
	void schrodigerSolverProb_instance1(unsigned num_qubits, std::complex<double>* H, unsigned dim, std::complex<double>* vectInitState, double time, std::complex<double>* result, std::complex<double>* Hadamard = NULL);

	template<typename T>
	void schrodigerSolverProb_instance2(unsigned num_qubits, std::complex<double>* H, unsigned dim, std::complex<double>* vectInitState, double time, std::complex<double>* result, std::complex<double>* Hadamard = NULL);

	template<typename T>
	vector<double> GetProbAmplitudes(std::complex<double>* vectState, unsigned number_of_qubits, unsigned atom_level, unsigned number_of_atom_level);

	template<typename T>
	void initZero(unsigned int N, std::complex<double> * A);

	template<typename T>
	void solver_diag_complex(unsigned int dim, std::complex<double> * matrix, std::complex<double> * result);

	template<typename T>
	void project_to_subspace(std::complex<double>* U, unsigned old_dim, std::complex<double>* U_proj, unsigned new_dim, unsigned new_elements_list[]);

	template<typename T>
	void solver_diag_real(unsigned int dim, double * matrix, std::complex<double> * result);

	template<typename T>
	void convert_hermitian_to_upper_triangle(std::complex<double> * a, unsigned dim);

	template<typename T>
	void print_comp_matrix(char* desc, int row, int col, std::complex<double> * a, int dim);

	template<typename T>
	void print_re_matrix(char* desc, int row, int col, double * a, int dim);

	template<typename T>
	double modulus(double x, double y);
}

//This will apply to both photon number n and transmon level j
template<typename T>
Operator<T> Utils::creation(const T& n) {
	std::deque<Projector<T>> creation_projs;

	for (T k = 0; k < n ; k++) {
		Projector<T> proj({k}, {k + 1}, sqrt(k+1));
		creation_projs.push_back(proj);
	}
	Operator<T> a_dagger(n + 1 , creation_projs);

	return a_dagger;
}

template<typename T>
Operator<T> Utils::annihilation(const T& n) {
	Operator<T> a = creation(n).dagger();

	return a;
}

template<typename T>
Operator<T> Utils::pauliI(const T& n) {
	std::deque<Projector<T>> I_projs;

	std::complex<double> weight(1, 0);
	for (T k = 0; k < n + 1; k++) {
		Projector<T> proj({ k }, { k }, weight);
		I_projs.push_back(proj);
	}

	Operator<T> I(n+1, I_projs);

	return I;
}

template<typename T>
Operator<T> Utils::pauliI() {
	std::deque<Projector<T>> I_projs;
	T k = 0;
	T j = 1;
	std::complex<double> weight(1, 0);

	Projector<T> proj1({ k }, { k }, weight);
	Projector<T> proj2({ j }, { j }, weight);

	I_projs.push_back(proj1);
	I_projs.push_back(proj2);

	Operator<T> I(2, I_projs);

	return I;
}

template<typename T>
Operator<T> Utils::pauliX() {
	std::deque<Projector<T>> sigmaX_projs;
	T k = 0;
	T j = 1;
	std::complex<double> weight(1, 0);

	Projector<T> proj1({ k }, { j }, weight);
	Projector<T> proj2({ j }, { k }, weight);

	sigmaX_projs.push_back(proj1);
	sigmaX_projs.push_back(proj2);

	Operator<T> X(2, sigmaX_projs);

	return X;
}

template<typename T>
Operator<T> Utils::pauliZ() {
	std::deque<Projector<T>> sigmaZ_projs;

	T k = 0;
	T j = 1;
	std::complex<double> weight1(1, 0);
	std::complex<double> weight2(-1, 0);

	Projector<T> proj1({ k }, { k }, weight1);
	Projector<T> proj2({ j }, { j }, weight2);

	sigmaZ_projs.push_back(proj1);
	sigmaZ_projs.push_back(proj2);

	Operator<T> Z(2, sigmaZ_projs);

	return Z;
}

template<typename T>
Operator<T> Utils::pauliY() {
	std::deque<Projector<T>> pauliY_projectors;

	T k = 0;
	T j = 1;
	std::complex<double> weight1(0, 1);
	std::complex<double> weight2(0, -1);

	Projector<T> proj1({ k }, { j }, weight1);
	Projector<T> proj2({ j }, { k }, weight2);

	pauliY_projectors.push_back(proj1);
	pauliY_projectors.push_back(proj2);

	Operator<T> pauliY(2, pauliY_projectors);

	return pauliY;
}

template<typename T>
Operator<T> Utils::pauliPlus() {
	std::deque<Projector<T>> projs;
	T k = 0;
	T j = 1;
	std::complex<double> weight(1, 0);

	Projector<T> proj1({ j }, { k }, weight);

	projs.push_back(proj1);

	Operator<T> opt(2, projs);

	return opt;
}

template<typename T>
Operator<T> Utils::pauliMinus() {
	std::deque<Projector<T>> projs;
	T k = 0;
	T j = 1;
	std::complex<double> weight(1, 0);

	Projector<T> proj1({ k }, { j }, weight);

	projs.push_back(proj1);

	Operator<T> opt(2, projs);

	return opt;
}

template<typename T>
Operator<T> Utils::Hadamard() {
	std::deque<Projector<T>> projs;
	T k = 0;
	T j = 1;
	std::complex<double> weight1(1/sqrt(2), 0);
	std::complex<double> weight2(-1/sqrt(2), 0);

	Projector<T> proj1({ k }, { k }, weight1);
	Projector<T> proj2({ j }, { k }, weight1);
	Projector<T> proj3({ k }, { j }, weight1);
	Projector<T> proj4({ j }, { j }, weight2);

	projs.push_back(proj1);
	projs.push_back(proj2);
	projs.push_back(proj3);
	projs.push_back(proj4);

	Operator<T> opt(2, projs);

	return opt;
}

template<typename T>
vector<double> Utils::GetProbAmplitudes(std::complex<double>* vectState, unsigned number_of_qubits, unsigned atom_level, unsigned number_of_atom_level) {
	
	vector<double> Q(number_of_qubits, 0.0);

	for (unsigned n = 0; n < number_of_qubits; n++)
	{
		for (unsigned i = atom_level * pow(number_of_atom_level, n); i < pow(number_of_atom_level, number_of_qubits);)
		{
			double sum = 0.0;

			for (unsigned j = i; j < (i + pow(number_of_atom_level, n)); j++)
			{
				std::complex<double> w = vectState[j];
				sum = sum + pow(abs(w), 2);
			}

			Q[n] = Q[n] + sum;
			i = i + pow(number_of_atom_level, n + 1);
		}
	}

	return Q;
}

template<typename T>
void Utils::initZero(unsigned int N, std::complex<double> * A)
{
	for (unsigned int n = 0; n < N; n++) {
		A[n]= std::complex<double>(0.0, 0.0);
	}

	return;
}

//Reference
//https://www.ibm.com/support/knowledgecenter/SSFHY8_5.2.0/com.ibm.cluster.essl.v5r2.essl100.doc/am5gr_hsspevx.htm
template<typename T>
void Utils::solver_diag_complex(unsigned int dim, std::complex<double> * matrix, std::complex<double> * result)
{
	// NOTE: the matrix is assumed to be symmetric, only the uppervalues are used!
	// Note: the matrix gets overwritten with the eigenvectors!
	unsigned int N = dim * dim;

	// 1) diagonalize the matrix: matrix = V * D * inv(V)
	//std::complex<double> * D = (std::complex<double> *)mkl_malloc(N * sizeof(std::complex<double>), 64);

	int matrix_layout = LAPACK_ROW_MAJOR;
	char jobz = 'V'; 
	char range = 'A'; 
	char uplo = 'U';
	lapack_int n = dim; 	
	
	/*printf("original hermitian matrix\n");
	for (unsigned i = 0; i < dim; i++) {
		for (unsigned j = 0; j < dim; j++) {
			std::complex<double> we = matrix[j + i * dim];
			printf("element: %u , weight: %lf + %lf i\n", (j + i * dim), we.real(), we.imag());
		}
		printf("\n");
	}*/
	Utils::convert_hermitian_to_upper_triangle<T>(matrix, dim);

	//printf("upper triangle matrix\n");
	//for (unsigned i = 0; i < dim; i++) {
	//	for (unsigned j = 0; j < dim; j++) {
	//		std::complex<double> we = matrix[j + i * dim];
	//		printf("element: %u , weight: %lf + %lf i\n", (j + i * dim), we.real(), we.imag());
	//	}
	//	printf("\n");
	//}

	//this a must be upper triangular (uplo = 'U') of input matrix, on exit the upper triangular including the diagonal of a is overwritten	
	MKL_Complex16 * a = reinterpret_cast <MKL_Complex16 *>(matrix);

	lapack_int lda = dim; 
	//double vl = 1; 
	//double vu = 100; 
	double vl = -100;
	double vu = 100;
	lapack_int il = 1; 
	lapack_int iu = 1; 
	double abstol = -1.0; 
	lapack_int m;
	lapack_int ldz = dim;

	double* w = (double *)mkl_malloc(dim * sizeof(double), 64);//Array, size at least max(1, n), contains the selected eigenvalues in ascending order, stored in w[0] to w[m - 1].
	lapack_complex_double* z = (lapack_complex_double *)mkl_malloc(N * sizeof(lapack_complex_double), 64);
	//output array If jobz = 'V', then if info = 0, the first m columns of z contain the orthonormal eigenvectors of the matrix A corresponding to the selected eigenvalues, 
	
	int* isuppz = (int *)mkl_malloc(dim * sizeof(int), 64);

	//with the i-th column of z holding the eigenvector associated with w[i - 1].
	
	//std::complex<double>* z = (std::complex<double> *)mkl_malloc(N * sizeof(std::complex<double>), 64);		
	//Array, size at least 2 *max(1, m).
	//The support of the eigenvectors in z, i.e., the indices indicating the nonzero elements in z.
	//The i - th eigenvector is nonzero only in elements isuppz[2i - 2] through isuppz[2i - 1].
	//Referenced only if eigenvectors are needed(jobz = 'V') and all eigenvalues are needed, that is, range = 'A' or range = 'I' and il = 1 and iu = n.
	
	int info = LAPACKE_zheevr(matrix_layout, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);

	//if info=0, successful
	if (info > 0) {
		cout << "The algorithm failed to compute eigenvalues" << endl;
		exit(1);
	}

	printf("w\n");
	for (unsigned i = 0; i< dim; i++) {	
		double we = w[i];
			printf("element: %u , weight: %lf \n", (i), we);
		printf("\n");
	}

	// 2) exponentiate -iD, elementwise, exp(-i*D)
	std::complex<double>* eig = (std::complex<double>*)mkl_malloc(N * sizeof(std::complex<double>), 64);

	initZero<T>(N, eig);
	for (unsigned int n = 0, i = 0; n < dim; n++, i += dim + 1)
	{
		eig[i] = std::complex<double>(cos(w[n]), -sin(w[n]));
	}

	// 3) multiply V * eig * transpose(V)
	std::complex<double>alpha(1.0, 0);
	std::complex<double>beta(0, 0);

	std::complex<double> * tmp = (std::complex<double> *)mkl_malloc(N * sizeof(std::complex<double>), 64);

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &alpha, z, dim, eig, dim, &beta, tmp, dim);

	printf("tmp\n");
	for (unsigned i = 0; i< dim; i++) {
		for (unsigned j = 0; j< dim; j++) {
			std::complex<double> w = tmp[j + i*dim];
			printf("element: %u , weight: %lf + %lf i\n", (j + i*dim), w.real(), w.imag());
		}
		printf("\n");
	}

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, dim, dim, dim, &alpha, tmp, dim, z, dim, &beta, result, dim);

	printf("result\n");
	for (unsigned i = 0; i< dim; i++) {
		for (unsigned j = 0; j< dim; j++) {
			std::complex<double> w = result[j + i*dim];
			printf("element: %u , weight: %lf + %lf i\n", (j + i*dim), w.real(), w.imag());
		}
		printf("\n");
	}

	mkl_free(w);
	mkl_free(eig);
	//mkl_free(a);
	mkl_free(tmp);
	mkl_free(z);
	mkl_free(isuppz);

	return;
}

template<typename T>
void Utils::schrodigerSolverProb_instance1(unsigned num_qubits, std::complex<double>* H, unsigned dim, std::complex<double>* vectInitState, double time, std::complex<double>* result, std::complex<double>* Hadamard)
{
	std::complex<double>alpha(1.0, 0);
	std::complex<double>beta(0, 0);

	unsigned computational_dim = pow(2, num_qubits);

	std::complex<double>* U_proj = (std::complex<double>*)mkl_malloc(computational_dim * computational_dim * sizeof(std::complex<double>), 64);
	std::complex<double>* U = (std::complex<double>*)mkl_malloc(dim * dim * sizeof(std::complex<double>), 64);
	std::complex<double>* R = (std::complex<double>*)mkl_malloc(dim * dim * sizeof(std::complex<double>), 64);

	for (unsigned i = 0; i < dim * dim; ++i) {
		R[i] = H[i] * time;
	}

	solver_diag_complex<T>(dim, R, U);

	if (Hadamard != NULL)
	{
		std::complex<double>* H_U_result = (std::complex<double>*)mkl_malloc(computational_dim * computational_dim * sizeof(std::complex<double>), 64);
		std::complex<double>* H_U_H_result = (std::complex<double>*)mkl_malloc(computational_dim * computational_dim * sizeof(std::complex<double>), 64);

		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			computational_dim, computational_dim, computational_dim, &alpha, Hadamard, computational_dim, U, computational_dim, &beta, H_U_result, computational_dim);

		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			computational_dim, computational_dim, computational_dim, &alpha, H_U_result, computational_dim, Hadamard, computational_dim, &beta, H_U_H_result, computational_dim);

		cblas_zcopy(computational_dim*computational_dim, H_U_H_result, 1, U, 1);
		mkl_free(H_U_result);
		mkl_free(H_U_H_result);
	}

	//https://software.intel.com/en-us/mkl-developer-reference-c-cblas-gemv
	cblas_zgemv(CblasRowMajor, CblasNoTrans, computational_dim, computational_dim, &alpha, U, computational_dim, vectInitState, 1, &beta, result, 1);


	mkl_free(R);
	mkl_free(U);
	mkl_free(U_proj);

	return;
}

//TODO: This method needs a lot of refactoring
template<typename T>
void Utils::schrodigerSolverProb_instance2(unsigned num_qubits, std::complex<double>* H, unsigned dim, std::complex<double>* vectInitState, double time, std::complex<double>* result, std::complex<double>* Hadamard)
{
	std::complex<double>alpha(1.0, 0);
	std::complex<double>beta(0, 0);

	//unsigned computational_dim = pow(2, num_qubits);
	//unsigned computational_space_list[] = { 0, 2, 6, 8 };

	unsigned computational_dim = 8;
	unsigned computational_space_list[] = { 0, 1, 2, 3, 6, 7, 8, 9};

	std::complex<double>* U_proj = (std::complex<double>*)mkl_malloc(computational_dim * computational_dim * sizeof(std::complex<double>), 64);
	std::complex<double>* U = (std::complex<double>*)mkl_malloc(dim * dim * sizeof(std::complex<double>), 64);
	std::complex<double>* R = (std::complex<double>*)mkl_malloc(dim * dim * sizeof(std::complex<double>), 64);

	for (unsigned i = 0; i < dim * dim; ++i) {
		R[i] = H[i] * time;
	}

	double* R_real = (double*)mkl_malloc(dim * dim * sizeof(double), 64);
	for (unsigned i = 0; i < dim * dim; i++) {
		R_real[i] = R[i].real();
	}
	solver_diag_real<T>(dim, R_real, U);

	Utils::project_to_subspace<T>(U, dim, U_proj, computational_dim, computational_space_list);

	//phase compensation
	std::complex<double>* V_temp = (std::complex<double>*)mkl_malloc(computational_dim * computational_dim * sizeof(std::complex<double>), 64);
	std::complex<double>* U_phase_compensated = (std::complex<double>*)mkl_malloc(computational_dim * computational_dim * sizeof(std::complex<double>), 64);

	std::complex<double> minus_i(0, -1);//-i	

	std::complex<double> beta0_0 = arg(U_proj[0]);
	std::complex<double> beta0_1 = arg(U_proj[1 + computational_dim]);
	std::complex<double> beta2_0 = arg(U_proj[2 + 2 * computational_dim]);
	std::complex<double> beta1_0 = arg(U_proj[4 + 4 * computational_dim]);

	std::complex<double> global_phase_compensation = exp(minus_i*beta0_0);

	Utils::initZero<unsigned>(computational_dim * computational_dim, V_temp);

	V_temp[0] = global_phase_compensation;

	V_temp[1 + 1 * computational_dim] = global_phase_compensation * exp(minus_i*beta0_1);

	V_temp[2 + 2 * computational_dim] = global_phase_compensation * exp(minus_i*beta2_0);
	V_temp[3 + 3 * computational_dim] = global_phase_compensation * exp(minus_i*(beta0_1 + beta2_0));

	V_temp[4 + 4 * computational_dim] = global_phase_compensation * exp(minus_i*beta1_0);
	V_temp[5 + 5 * computational_dim] = global_phase_compensation * exp(minus_i*(beta1_0 + beta0_1));

	V_temp[6 + 6 * computational_dim] = global_phase_compensation * exp(minus_i*(beta1_0 + beta2_0));
	V_temp[7 + 7 * computational_dim] = global_phase_compensation * exp(minus_i*(beta1_0 + beta2_0 + beta0_1));

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, computational_dim, computational_dim, computational_dim, &alpha, U_proj, computational_dim, V_temp, computational_dim, &beta, U_phase_compensated, computational_dim);

	printf("U_phase_compensated matrix\n");
	for (unsigned i = 0; i < computational_dim; i++) {
		for (unsigned j = 0; j < computational_dim; j++) {
			std::complex<double> w = U_phase_compensated[j + i * computational_dim];
			printf("element: %u , weight: %lf + %lf i\n", (j + i * computational_dim), w.real(), w.imag());
		}
		printf("\n");
	}

	if (Hadamard != NULL)
	{
		std::complex<double>* H_U_result = (std::complex<double>*)mkl_malloc(computational_dim * computational_dim * sizeof(std::complex<double>), 64);
		std::complex<double>* H_U_H_result = (std::complex<double>*)mkl_malloc(computational_dim * computational_dim * sizeof(std::complex<double>), 64);

		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			computational_dim, computational_dim, computational_dim, &alpha, Hadamard, computational_dim, U_phase_compensated, computational_dim, &beta, H_U_result, computational_dim);

		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			computational_dim, computational_dim, computational_dim, &alpha, H_U_result, computational_dim, Hadamard, computational_dim, &beta, H_U_H_result, computational_dim);

		cblas_zcopy(computational_dim*computational_dim, H_U_H_result, 1, U_phase_compensated, 1);
		mkl_free(H_U_result);
		mkl_free(H_U_H_result);
	}

	//https://software.intel.com/en-us/mkl-developer-reference-c-cblas-gemv
	cblas_zgemv(CblasRowMajor, CblasNoTrans, computational_dim, computational_dim, &alpha, U_phase_compensated, computational_dim, vectInitState, 1, &beta, result, 1);

	mkl_free(R);
	mkl_free(U);
	mkl_free(U_proj);
	mkl_free(V_temp);
	mkl_free(U_phase_compensated);

	mkl_free(R_real);

	return;
}

template<typename T>
void Utils::project_to_subspace(std::complex<double>* U,unsigned old_dim, std::complex<double>* U_proj, unsigned new_dim, unsigned new_elements_list[])
{
	 for (int i = 0; i < new_dim; i++)
	 {
		 for (int j = 0; j < new_dim; j++)
		 {
			 U_proj[j + i * new_dim] = U[new_elements_list[j] + new_elements_list[i] * old_dim];
		 }
	 }
}

template<typename T>
void Utils::solver_diag_real(unsigned int dim, double * matrix, std::complex<double> * result)
{
	// NOTE: the matrix is assumed to be real symmetric, only the uppervalues are used!
	// Note: the matrix gets overwritten with the eigenvectors!
	unsigned int N = dim * dim;

	// 1) diagonalize the matrix: matrix = V * D * inv(V)
	double * D = (double*)mkl_malloc(dim * sizeof(double), 64);

	int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', dim, matrix, dim, D);

	if (info != 0) {
		cout << "The algorithm failed to compute eigenvalues" << endl;
		exit(1);
	}

	// 2) exponentiate -iD, elementwise, exp(-i*D)
	std::complex<double>* eig = (std::complex<double>*)mkl_malloc(N * sizeof(std::complex<double>), 64);

	initZero<T>(N, eig);
	for (unsigned int n = 0, i = 0; n < dim; n++, i += dim + 1)
	{
		eig[i] = std::complex<double>(cos(D[n]), -sin(D[n]));
	}

	// 3) multiply V * eig * inv(V)
	std::complex<double> zero = { 0.0, 0.0 };
	std::complex<double> one = { 1.0, 0.0 };

	std::complex<double>* tmp = (std::complex<double>*)mkl_malloc(N * sizeof(std::complex<double>), 64);
	std::complex<double>* V = (std::complex<double>*)mkl_malloc(N * sizeof(std::complex<double>), 64);

	for (unsigned int n = 0; n < N; n++) 
	{
		V[n] = std::complex<double>(matrix[n], 0.0);
	}

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &one, V, dim, eig, dim, &zero, tmp, dim);
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, dim, dim, dim, &one, tmp, dim, V, dim, &zero, result, dim);

	mkl_free(D);
	mkl_free(eig);
	mkl_free(tmp);
	mkl_free(V);
}

template<typename T>
void Utils::convert_hermitian_to_upper_triangle(std::complex<double> * a, unsigned dim)
{
	std::complex<double>zero(0, 0);
	for (int i = 1; i < dim; i++) 
	{
		for (int j = 0; j < i; j++)
		{
			a[j + i * dim] = zero;
		}
	}
	return;
}

template<typename T>
void Utils::print_comp_matrix(char* desc, int row, int col, std::complex<double> * a, int dim)
{
	printf("\n %s\n", desc);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++)
			printf(" (%6.2f,%6.2f)", a[i + j * dim].real(), a[i + j * dim].imag());
		printf("\n");
	}
}

template<typename T>
void Utils::print_re_matrix(char* desc, int row, int col, double * a, int dim) 
{
	int i, j;
	printf("\n %s\n", desc);
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) printf(" %6.2f", a[i + j * dim]);
		printf("\n");
	}
}

template<typename T>
double Utils::modulus(double x, double y) 
{
	double result = x;
	if (x > y || x < -y)
	{
		double m = x / y;
		double r = m * y;
		result = x - r;
	}

	return result;
}

#endif