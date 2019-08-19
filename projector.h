#pragma once

#ifndef __PROJECTOR_H
#define __PROJECTOR_H

#include "dirac.h"
#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

template <typename T> class Projector {
public:
	Bra<T> bra;
	Ket<T> ket;
	std::complex<double> weight;

	// Default Constructor                                                                                                                                                      
	Projector();

	Projector(Bra<T> bra, Ket<T> ket, std::complex<double> weight);

	Projector(std::vector<T> bra_elements, std::vector<T> ket_elements, std::complex<double> weight);

	// Copy Constructor    
	Projector(const Projector<T>& rhs);
	virtual ~Projector();

	// Operator overloading, for "standard" mathematical dirac operations                                                                                                                                                          
	Projector<T>& operator=(const Projector<T>& rhs);

	// Dirac mathematical operations                                                                                                                                                                                               
	std::deque<Projector<T>> operator+(const Projector<T>& rhs);
	std::deque<Projector<T>> operator-(const Projector<T>& rhs);
	Projector<T> operator*(const Projector<T>& rhs);

	Projector<T> inner(const Projector<T>& rhs);
	//Projector<T> innerleft(const Projector<T>& lhs);
	Projector<T> transpose();
	Projector<T> conjugate();

	// Bra/scalar multiplication
	Projector<T> operator*(const T& rhs);
	// Bra/scalar division                                                                                                                                                     
	Projector<T> operator/(const T& rhs);

	Ket<T> applyTo(const Ket<T>& rhs);

	// Ket/vector operations                                                                                                                                                                                                     
	//std::vector<T> operator*(const std::vector<T>& rhs);
	//std::vector<T> diag_vec();
	bool matchKet(const Projector<T> rhs);
	bool matchBra(const Projector<T> rhs);
	bool match(const Projector<T> rhs);
	void print();
};

#include "projector.cpp"

#endif