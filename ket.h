#pragma once

#ifndef __KET_H
#define __KET_H

#include "dirac.h"
#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

template <typename T> class Bra;
template <typename T> class Projector;
template <typename T> class State;

template <typename T> class Ket {
public:
	std::vector<T> elements;
	std::complex<double> weight;
	unsigned length;

	Ket();
	Ket(unsigned length, const std::vector<T> elements, std::complex<double> weight);
	Ket(const Ket<T>& rhs);
	virtual ~Ket();

	// Operator overloading, for "standard" mathematical dirac operations                                                                                                                                                          
	Ket<T>& operator=(const Ket<T>& rhs);

	bool operator==(const Ket<T>& rhs);

	// Dirac mathematical operations                                                                                                                                                                                               
	State<T> operator+(const Ket<T>& rhs);
	std::deque<Ket<T>> operator-(const Ket<T>& rhs);
	// Left tensor product of this ket and another ket
	Ket<T> operator*(const Ket<T>& rhs);

	// Inner product of this Ket with a bra applied from left
	std::complex<double> inner(const Bra<T>& lhs);
	Projector<T> outer(const Bra<T>& rhs);
	Bra<T> transpose();
	Ket<T> conjugate();

	// Bra/scalar multiplication
	Ket<T> operator*(const T& rhs);
	// Bra/scalar division                                                                                                                                                     
	Ket<T> operator/(const T& rhs);

	// Ket/vector operations                                                                                                                                                                                                     
	//std::vector<T> operator*(const std::vector<T>& rhs);
	//std::vector<T> diag_vec();

	// Access the individual element                                                                                                                                                                                               
	T& operator()(const unsigned& element);
	const T& operator()(const unsigned& element) const;

	// Access the elements size                                                                                                                                                                                              
	unsigned get_elements() const;

	std::pair<unsigned, std::complex<double>> asDecimal(unsigned base);

	void print();
};

#include "ket.cpp"

#endif
