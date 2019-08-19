#pragma once

#ifndef __STATE_H
#define __STATE_H

#include "dirac.h"
#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

template <typename T> class State {
public:
	std::deque<Ket<T>> elements;
	unsigned size;

	State();
	State(unsigned size, const std::deque<Ket<T>> elements);
	State(const State<T>& rhs);
	virtual ~State();

	// Operator overloading, for "standard" mathematical dirac operations                                                                                                                                                          
	State<T>& operator=(const State<T>& rhs);

	State<T> operator+(const Ket<T>& rhs);

	// State/scalar multiplication
	State<T> operator*(const T& rhs);
	// State/scalar division                                                                                                                                                     
	State<T> operator/(const T& rhs);

	// Access the individual element                                                                                                                                                                                               
	std::deque<Ket<T>>& operator()(const unsigned& element);
	const std::deque<Ket<T>>& operator()(const unsigned& element) const;

	// Get the number of elements of operator                                                                                                                                                                                           
	unsigned get_elements() const;

	std::complex<double>* asVector(unsigned base);

	void print();
};

#include "state.cpp"

#endif
