#ifndef __BRA_CPP
#define __BRA_CPP

#include "bra.h"

#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

#pragma region Bra

// Default Constructor                                                                                                                                                      
template<typename T>
Bra<T>::Bra() {}

// Parameter Constructor                                                                                                                                                      
template<typename T>
Bra<T>::Bra(unsigned _length, const std::vector<T> _elements, std::complex<double> _weight) {
	length = _length;
	elements.resize(_length);
	elements = _elements;
	weight = _weight;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
Bra<T>::Bra(const Bra<T>& rhs) {
	length = rhs.get_elements();
	elements = rhs.elements;
	weight = rhs.weight;
}

// (Virtual) Destructor                                                                                                                                                       
template<typename T>
Bra<T>::~Bra() {}

// Assignment Operator                                                                                                                                                        
template<typename T>
Bra<T>& Bra<T>::operator=(const Bra<T>& rhs) {
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

// Addition of two Bras                                                                                                                                                 
template<typename T>
Bra<T>& Bra<T>::operator+(const Bra<T>& rhs) {

	for (unsigned i = 0; i<length; i++) {
		//if two bras have the same values, then the weights can be added
		if (elements[i] != rhs(i)) {
			string message = "Two Bras with different values cannot be added";
			printf(message);
			throw invalid_argument(message);
		}
	}

	this->weight += rhs.weight;

	return *this;
}

// Subtraction of this Bra and another                                                                                                                                     
template<typename T>
Bra<T>& Bra<T>::operator-(const Bra<T>& rhs) {

	for (unsigned i = 0; i<length; i++) {
		//if two bras have the same values, then the weights can be added
		if (elements[i] != rhs(i)) {
			string message = "Two Bras with different values cannot be subtracted";
			printf(message);
			throw invalid_argument(message);
		}
	}

	this->weight -= rhs.weight;

	return *this;
}

// Left tensor product of this Bra and another                                                                                                                              
template<typename T>
Bra<T> Bra<T>::operator*(const Bra<T>& rhs) {
	std::vector<T> new_elements = this->elements;

	for (unsigned i = 0; i < rhs.get_elements(); i++) {
		new_elements.push_back(rhs(i));
	}
	std::complex<double> new_weight = this->weight * rhs.weight;
	unsigned new_length = this->length + rhs.length;
	Bra result(new_length, new_elements, new_weight);

	return result;
}

// Inner product of this Bra with a ket applied to right
//if the bra value is equal to ket value, result is multiplication of weights, else 0
template<typename T>
std::complex<double> Bra<T>::inner(const Ket<T>& rhs) {
	for (unsigned i = 0; i<length; i++) {
		if (elements[i] != rhs(i)) {
			std::complex<double> zero;
			return zero;
		}
	}
	std::complex<double> result = this->weight * rhs.weight;

	return result;
}

// Return the transpose of this bra                                                                                                                                       
template<typename T>
Ket<T> Bra<T>::transpose() {
	Ket<T> result(this->length, this->elements, this->weight);

	return result;
}

// Return the complex conjugate of this bra                                                                                                                                       
template<typename T>
Bra<T> Bra<T>::conjugate() {
	this->weight = std::conj(this->weight);

	return *this;
}

// Outer product of this Bra with a ket applied from the left	
#pragma warning //TODO: check if the length of the ket and bra should be equal!
template<typename T>
Projector<T> Bra<T>::outer(const Ket<T>& lhs) {
	std::vector<T> ket_elements = lhs.elements;
	std::vector<T> bra_elements = this->elements;

	Ket<T> new_ket(this->length, ket_elements, 1);
	Bra<T> new_bra(this->length, bra_elements, 1);
	std::complex<double> new_weight = this->weight * rhs.weight;

	Projector<T> new_projector(new_bra, new_ket, new_weight);

	return new_projector;
}


// Bra/scalar multiplication                                                                                                                                               
template<typename T>
Bra<T> Bra<T>::operator*(const T& rhs) {
	this->weight *= rhs

		return *this;
}

// Bra/scalar division                                                                                                                                                     
template<typename T>
Bra<T> Bra<T>::operator/(const T& rhs) {
	this->weight /= rhs

		return *this;
}

// Multiply a bra with a projector                                                                                                                                            
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
T& Bra<T>::operator()(const unsigned& index) {
	return this->elements[index];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& Bra<T>::operator()(const unsigned& index) const {
	return this->elements[index];
}

// Get the number of elements of Bra                                                                                                                                       
template<typename T>
unsigned Bra<T>::get_elements() const {
	return this->elements.size();
}

template<typename T>
std::pair<unsigned, std::complex<double>> Bra<T>::asDecimal(unsigned base) {
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
void Bra<T>::print() {
	printf("Bra weight: %lf + %lf i\n", this->weight.real(), this->weight.imag());
	printf("Bra elements:");
	for (int i = 0; i < this->get_elements(); i++) {
		printf("%u", this->elements[i]);
	}
	printf("\n");
}
#pragma endregion

#endif