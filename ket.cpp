#ifndef __KET_CPP
#define __KET_CPP

#include "ket.h"

#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

#pragma region Ket

// Default Constructor                                                                                                                                                      
template<typename T>
Ket<T>::Ket() {}

// Parameter Constructor                                                                                                                                                      
template<typename T>
Ket<T>::Ket(unsigned _length, std::vector<T> _elements, std::complex<double> _weight) {
	length = _length;
	elements.resize(_length);
	elements = _elements;
	weight = _weight;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
Ket<T>::Ket(const Ket<T>& rhs) {
	length = rhs.get_elements();
	elements = rhs.elements;
	weight = rhs.weight;
}

// (Virtual) Destructor                                                                                                                                                       
template<typename T>
Ket<T>::~Ket() {}

// Assignment Operator                                                                                                                                                        
template<typename T>
Ket<T>& Ket<T>::operator=(const Ket<T>& rhs) {
	if (&rhs == this)
		return *this;

	unsigned rhs_lenght = rhs.get_elements();

	elements.resize(rhs_lenght);

	for (unsigned i = 0; i < rhs_lenght; i++) {
		elements[i] = rhs(i);
	}

	length = rhs.length;
	weight = rhs.weight;

	return *this;
}


// Equality Operator                                                                                                                                                        
template<typename T>
bool Ket<T>::operator==(const Ket<T>& rhs) {
	if (&rhs == this)
		return true;
	else
		return false;
	//unsigned rhs_lenght = rhs.get_elements();

	//elements.resize(rhs_lenght);

	//for (unsigned i = 0; i < rhs_lenght; i++) {
	//	elements[i] = rhs(i);
	//}

	//length = rhs.length;
	//weight = rhs.weight;

	//return *this;
}
// Addition of this ket with another ket                                                                                                                                                 
template<typename T>
State<T> Ket<T>::operator+(const Ket<T>& rhs) {
	std::deque<Ket<T>> result;
	bool equalBase = true;
	for (unsigned i = 0; i<length; i++) {
		if (elements[i] != rhs(i)) {
			equalBase = false;
			break;
		}
	}
	if (equalBase) {
		this->weight += rhs.weight;		
	}
	else
	{
		result.push_back(rhs);
	}
	result.push_back(this);

	State<T> result_state(result.size(), result);
	return result_state;
}
//TODO: fix this
// Subtraction of this Ket and another                                                                                                                                     
template<typename T>
std::deque<Ket<T>> Ket<T>::operator-(const Ket<T>& rhs) {
	std::deque<Ket<T>> result;
	bool equalBase = true;
	for (unsigned i = 0; i<length; i++) {
		if (elements[i] != rhs(i)) {
			equalBase = false;
			break;
		}
	}
	if (equalBase) {
		this->weight -= rhs.weight;
	}
	else
	{
		result.push_back(rhs);
	}
	result.push_back(this);

	return result;
}

// Left tensor product of this ket and another ket                                                                                                                             
template<typename T>
Ket<T> Ket<T>::operator*(const Ket<T>& rhs) {
	std::vector<T> new_elements = this->elements;

	for (unsigned i = 0; i < rhs.get_elements(); i++) {
		new_elements.push_back(rhs(i));
	}
	std::complex<double> new_weight = this->weight * rhs.weight;
	unsigned new_length = this->length + rhs.length;
	Ket<T> result(new_length, new_elements, new_weight);

	return result;
}

// Inner product of this Ket with a bra applied from the left
template<typename T>
std::complex<double> Ket<T>::inner(const Bra<T>& lhs) {
	for (unsigned i = 0; i<length; i++) {
		//if the ket value (elements) is equal to bra value (elements), result is multiplication of weights, else 0
		if (elements[i] != rhs(i))
		{
			std::complex<double> zero;
			return zero;
		}
	}
	double result = this->weight * rhs.weight;

	return result;
}

//// Outer product of this Ket with a bra applied to the right	
////check the length of the ket and bra be equal
template<typename T>
Projector<T> Ket<T>::outer(const Bra<T>& rhs) {
	std::vector<T> bra_elements = rhs.elements;
	std::vector<T> ket_elements = this->elements;

	Ket<T> new_ket(this->length, ket_elements, 1);
	Bra<T> new_bra(this->length, bra_elements, 1);
	std::complex<double> new_weight = this->weight * rhs.weight;

	Projector<T> new_projector(new_bra, new_ket, new_weight);

	return new_projector;
}

// Return the transpose of this Ket                                                                                                                                       
template<typename T>
Bra<T> Ket<T>::transpose() {
	Bra<T> result(this->length, this->elements, this->weight);

	return result;
}

// Return the complex conjugate of this bra                                                                                                                                       
template<typename T>
Ket<T> Ket<T>::conjugate() {
	this->weight = std::conj(this->weight);

	return *this;
}

// Bra/scalar multiplication                                                                                                                                               
template<typename T>
Ket<T> Ket<T>::operator*(const T& rhs) {
	this->weight *= rhs

		return *this;
}

// Bra/scalar division                                                                                                                                                     
template<typename T>
Ket<T> Ket<T>::operator/(const T& rhs) {
	this->weight /= rhs

		return *this;
}

// Multiply a Ket with a projector                                                                                                                                            
//template<typename T>
//std::vector<T> QSMatrix<T>::operator*(const std::vector<T>& rhs) {
//	std::vector<T> result(rhs.size(), 0.0);
//
//	for (unsigned i = 0; i<rows; i++) {
//		for (unsigned j = 0; j<cols; j++) {
//			result[i] += this->mat[i][j] * rhs[j];
//		}
//	}
//
//	return result;
//}

// Access the individual elements                                                                                                                                             
template<typename T>
T& Ket<T>::operator()(const unsigned& index) {

	return this->elements[index];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& Ket<T>::operator()(const unsigned& index) const {

	return this->elements[index];
}

// Get the number of elements of Bra                                                                                                                                       
template<typename T>
unsigned Ket<T>::get_elements() const {
	return this->elements.size();
}

template<typename T>
std::pair<unsigned, std::complex<double>> Ket<T>::asDecimal(unsigned base) {
	unsigned result = 0;
	int j = 0;
	unsigned element_size = this->elements.size();

	for (int i = element_size - 1; i >= 0; i--) {

		unsigned value = this->elements[j];
		result = value * pow(base, i) + result;
		j++;
	}

	return std::pair<unsigned, std::complex<double>>(result, this->weight);
}

template<typename T>
void Ket<T>::print() {
	printf("Ket weight: %lf + %lf i\n", this->weight.real(), this->weight.imag());
	printf("Ket elements:");
	for (int i = 0; i < this->get_elements(); i++) {
		printf("%u", this->elements[i]);
	}
	printf("\n");
}
#pragma endregion

#endif