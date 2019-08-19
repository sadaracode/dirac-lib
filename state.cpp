#ifndef __STATE_CPP
#define __STATE_CPP

#include "state.h"

#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

#pragma region State

// Default Constructor                                                                                                                                                      
template<typename T>
State<T>::State() {}

// Parameter Constructor                                                                                                                                                      
template<typename T>
State<T>::State(unsigned _size, const std::deque<Ket<T>> _elements) {
	elements.resize(_elements.size());
	elements = _elements;
	size = _size;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
State<T>::State(const State<T>& rhs) {
	elements = rhs.elements;
	size = rhs.size;
}

// (Virtual) Destructor                                                                                                                                                       
template<typename T>
State<T>::~State() {}

// Assignment Operator                                                                                                                                                        
template<typename T>
State<T>& State<T>::operator=(const State<T>& rhs) {
	if (&rhs == this)
		return *this;

	this->elements = rhs.elements;
	this->size = rhs.size;
	return *this;
}

// Addition of this state with a ket                                                                                                                                                 
template<typename T>
State<T> State<T>::operator+(const Ket<T>& rhs) {
	std::complex<double> w(1);
	Ket<T> rhs_basis(rhs.length, rhs.elements, w);

	std::deque<Ket<T>> this_kets = this->elements;
	bool found = false;

	for (unsigned i = 0; i < this->get_elements(); i++)
	{
		Ket<T> ket = this_kets.front();
		if (ket.elements == rhs_basis.elements) 
		{
			ket.weight += rhs.weight;
			found = true;
			break;
		}
		this_kets.pop_front();
		this_kets.push_back(ket);
	}

	if (!found) {
		this->elements.push_back(rhs);
	}

	return *this;
}

// Operator/scalar multiplication                                                                                                                                               
template<typename T>
State<T> State<T>::operator*(const T& rhs) {
	for (unsigned i = 0; i < this->get_elements(); i++) {
		this->elements(i).weight *= rhs
	}

	return *this;
}

// Operator/scalar division                                                                                                                                                     
template<typename T>
State<T> State<T>::operator/(const T& rhs) {
	for (unsigned i = 0; i < this->get_elements(); i++) {
		this->elements(i).weight /= rhs
	}

	return *this;
}

// Access the individual elements                                                                                                                                             
template<typename T>
std::deque<Ket<T>>& State<T>::operator()(const unsigned& index) {

	return this->elements[index];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const std::deque<Ket<T>>& State<T>::operator()(const unsigned& index) const {

	return this->elements[index];
}

// Get the number of elements of Operator                                                                                                                                       
template<typename T>
unsigned State<T>::get_elements() const {

	return this->elements.size();
}

template<typename T>
std::complex<double>* State<T>::asVector(unsigned base) {

	std::complex<double> zero;
	//std::vector<std::complex<double>> result(this->size, zero);
	unsigned dim = this->size;
	std::complex<double>* result = (std::complex<double>*)mkl_malloc(dim * sizeof(std::complex<double>), 64);

	for (unsigned i = 0; i < dim; i++) {
		result[i] = zero;
	}

	std::deque<Ket<T>>::const_iterator it = this->elements.begin();
	while (it != this->elements.end()) {
		//for test
		Ket<T> k = *it;
		std::pair<unsigned, std::complex<double>> decimalKet = k.asDecimal(base);
		result[decimalKet.first] = decimalKet.second;
		++it;
	}
	/*for (unsigned i = 0; i < this->get_elements(); i++) {
		Ket<T> ket = this->elements(i);
		std::pair<unsigned, std::complex<double>> decimalKet = ket.asDecimal(base);
		result[decimalKet.first] = decimalKet.second;
	}*/

	return result;
}

template<typename T>
void State<T>::print() {
	printf("Start printing of this State \n\n");
	for (unsigned i = 0; i < this->get_elements(); i++) {
		//TODO: remove variable after test
		Ket<T> k = this->elements[i];
		k.print();
	}
	printf("End printing of this State \n\n");
}
#pragma endregion

#endif