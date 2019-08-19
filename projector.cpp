#ifndef __PROJECTOR_CPP
#define __PROJECTOR_CPP

#include "projector.h"

#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

#pragma region Projector

// Default Constructor                                                                                                                                                      
template<typename T>
Projector<T>::Projector() {}

// Parameter Constructor 
template<typename T>
Projector<T>::Projector(Bra<T> _bra, Ket<T> _ket, std::complex<double> _weight) {
	bra = _bra;
	ket = _ket;
	weight = _weight;
}

template<typename T>
Projector<T>::Projector(std::vector<T> _bra_elements, std::vector<T> _ket_elements, std::complex<double> _weight) {
	Bra<T> _bra(1, _bra_elements, 1.0);
	Ket<T> _ket(1, _ket_elements, 1.0);
	bra = _bra;
	ket = _ket;
	weight = _weight;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
Projector<T>::Projector(const Projector<T>& rhs) {
	bra = rhs.bra;
	ket = rhs.ket;
	weight = rhs.weight;
}


// (Virtual) Destructor                                                                                                                                                       
template<typename T>
Projector<T>::~Projector() {}

// Assignment Operator                                                                                                                                                        
template<typename T>
Projector<T>& Projector<T>::operator=(const Projector<T>& rhs) {
	if (&rhs == this)
		return *this;

	this->bra = rhs.bra;
	this->ket = rhs.ket;
	this->weight = rhs.weight;

	return *this;
}

// Addition of two Projectors                                                                                                                                                 
template<typename T>
std::deque<Projector<T>> Projector<T>::operator+(const Projector<T>& rhs) {
	std::deque<Projector<T>> result;

	if (this->match(rhs)) {
		this->weight += rhs.weight;
	}
	else
	{
		result.push_back(rhs);
	}

	result.push_back(this);

	return result;
}

#pragma warning: TODO: update subtraction

// Subtraction of this Projector and another                                                                                                                                     
template<typename T>
std::deque<Projector<T>> Projector<T>::operator-(const Projector<T>& rhs) {
	std::deque<Projector<T>> result;
	if (this->match(rhs)) {
		this->weight -= rhs.weight;
	}
	else
	{
		result.push_back(rhs);
	}

	result.push_back(this);

	return result;
}

// Left tensor product of this Projector and another                                                                                                                              
template<typename T>
Projector<T> Projector<T>::operator*(const Projector<T>& rhs) {
	Bra<T> result_bra = this->bra * rhs.bra;
	Ket<T> result_ket = this->ket * rhs.ket;
	std::complex<double> result_weight = this->weight * rhs.weight;

	Projector<T> result(result_bra, result_ket, result_weight);

	return result;
}

// Inner product of this Projector with another Projector applied to right	
// If bra_ket_result is zero, the projector weight will be zero.
template<typename T>
Projector<T> Projector<T>::inner(const Projector<T>& rhs) {
	std::complex<double> bra_ket_result = this->bra.inner(rhs.ket);
	std::complex<double> result_weight = bra_ket_result * this->weight * rhs.weight;

	Projector<T> result(rhs.bra, this->ket, result_weight);

	return result;
}

// Return the transpose of this projector                                                                                                                                       
template<typename T>
Projector<T> Projector<T>::transpose() {
	Projector<T> result(this->ket.elements, this->bra.elements, this->weight);

	return result;
}

// Return the conjugate of this projector                                                                                                                                       
template<typename T>
Projector<T> Projector<T>::conjugate() {
	this->weight = std::conj(this->weight);

	return *this;
}

//TODO: implement dagger

// Projector/scalar multiplication                                                                                                                                               
template<typename T>
Projector<T> Projector<T>::operator*(const T& rhs) {
	this->weight *= rhs

		return *this;
}

// Projector/scalar division                                                                                                                                                     
template<typename T>
Projector<T> Projector<T>::operator/(const T& rhs) {
	this->weight /= rhs

		return *this;
}

// Apply a Projector to a ket could result weight = 0                                                                                                                                            
template<typename T>
Ket<T> Projector<T>::applyTo(const Ket<T>& rhs) {
	std::complex<double> result_weight = this->bra.inner(rhs) * this->weight;
	std::vector<T> elements = this->ket.elements;
	Ket<T> result(rhs.length, elements, result_weight);

	return result;
}

// Returns true if this projecr ket basis matches the rhs's ket basis                                                                                                                                          
template<typename T>
bool Projector<T>::matchKet(const Projector<T> rhs) {

	return this->ket.elements == rhs.ket.elements;
}

// Returns true if this projecr bra basis matches the rhs's bra basis                                                                                                                                          
template<typename T>
bool Projector<T>::matchBra(const Projector<T> rhs) {

	return this->bra.elements == rhs.bra.elements;
}

// Returns true if this projecr bra basis matches the rhs's bra basis                                                                                                                                          
template<typename T>
bool Projector<T>::match(const Projector<T> rhs) {

	return matchKet(rhs) && matchBra(rhs);
}

template<typename T>
void Projector<T>::print() {
	printf("Projector weight: %lf + %lf i\n", this->weight.real(), this->weight.imag()); 

	this->ket.print();
	printf("\n");

	this->bra.print();
	printf("\n");
}
#pragma endregion

#endif
