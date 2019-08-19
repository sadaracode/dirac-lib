#pragma once

#ifndef __BRA_H
#define __BRA_H

#include  "dirac.h"
#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

template <typename T> class Ket ;
template <typename T> class Projector ;

template <typename T> class Bra {
public:
	std::vector<T> elements;
	std::complex<double> weight;
	unsigned length;

	Bra();
	Bra(unsigned length, const std::vector<T> elements, std::complex<double> weight);
	Bra(const Bra<T>& rhs);
	virtual ~Bra();

	// Operator overloading, for "standard" mathematical dirac operations                                                                                                                                                          
	Bra<T>& operator=(const Bra<T>& rhs);

	// Dirac mathematical operations                                                                                                                                                                                               
	Bra<T>& operator+(const Bra<T>& rhs);
	Bra<T>& operator-(const Bra<T>& rhs);
	Bra<T> operator*(const Bra<T>& rhs);

	// Inner product of this Bra with a ket
	std::complex<double> inner(const Ket<T>& rhs);
	//unsigned inner(const std::vector<T>& ket_elements, unsigned ket_weight);
	Ket<T> transpose();
	Bra<T> conjugate();

	Projector<T> outer(const Ket<T>& lhs);
	//Projector<T&> outer(const std::vector<T&> ket_elements, unsigned ket_weight, T&);

	// Bra/scalar multiplication
	Bra<T> operator*(const T& rhs);
	// Bra/scalar division                                                                                                                                                     
	Bra<T> operator/(const T& rhs);

	// Bra/vector operations                                                                                                                                                                                                     
	//std::vector<T> operator*(const std::vector<T>& rhs);
	//std::vector<T> diag_vec();

	// Access the individual element                                                                                                                                                                                               
	T& operator()(const unsigned& element);
	const T& operator()(const unsigned& element) const;

	// Access the row and column sizes                                                                                                                                                                                              
	unsigned get_elements() const;

	std::pair<unsigned, std::complex<double>> asDecimal(unsigned base);

	void print();
};

#include "bra.cpp"

#endif