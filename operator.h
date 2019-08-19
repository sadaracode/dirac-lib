#pragma once

#ifndef __OPERATOR_H
#define __OPERATOR_H

#include "dirac.h"

#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>
#define mkl_malloc MKL_malloc

template <typename T> class State;

template <typename T> class Operator {
public:
	std::deque<Projector<T>> elements;
	unsigned size; //If the Operator is like a n by n matrix, the size is n

	Operator();
	Operator(unsigned size, const std::deque<Projector<T>> elements);
	Operator(const Operator<T>& rhs);
	virtual ~Operator();

	// Operator overloading, for "standard" mathematical dirac operations                                                                                                                                                          
	Operator<T>& operator=(const Operator<T>& rhs);

	// Dirac mathematical operations                                                                                                                                                                                               
	Operator<T> operator+(const Operator<T>& rhs);
	Operator<T> operator-(const Operator<T>& rhs);
	Operator<T> operator*(const Operator<T>& rhs);

	// Inner product of this Operator with a Operator
	Operator<T> inner(const Operator<T>& rhs);
	const Operator<T> transpose();
	const Operator<T> conjugate();
	const Operator<T> dagger();

	//TODO: Commutator

	//Projector<T> outer(const Ket<T>& lhs);
	//Projector<T&> outer(const std::vector<T&> ket_elements, unsigned ket_weight, T&);

	// Bra/scalar multiplication
	Operator<T> operator*(const std::complex<double> rhs);
	// Bra/scalar division                                                                                                                                                     
	Operator<T> operator/(const std::complex<double> rhs);

	// Bra/vector operations                                                                                                                                                                                                     
	//std::vector<T> operator*(const std::vector<T>& rhs);
	//std::vector<T> diag_vec();

	// Access the individual element                                                                                                                                                                                               
	T& operator()(const unsigned& element);
	const T& operator()(const unsigned& element) const;

	// Get the number of elements of operator                                                                                                                                                                                           
	unsigned get_elements() const;

	State<T> applyTo(const State<T>& rhs);

	void asMatrix(unsigned num_qubits, unsigned base, std::complex<double> *R);

	double * asRealMatrix(unsigned base);

	void print();
};

#include "operator.cpp"

#endif