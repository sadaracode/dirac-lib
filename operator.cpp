#ifndef __OPERATOR_CPP
#define __OPERATOR_CPP

#include "operator.h"

#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

// May change data structure of Operaqtor to map of key(projector.bra), value (projector)
#pragma region Operator

// Default Constructor                                                                                                                                                      
template<typename T>
Operator<T>::Operator() {}

// Parameter Constructor                                                                                                                                                      
template<typename T>
Operator<T>::Operator(unsigned _size, const std::deque<Projector<T>> _elements) {
	elements.resize(_elements.size());
	elements = _elements;
	size = _size;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
Operator<T>::Operator(const Operator<T>& rhs) {
	elements = rhs.elements;
	size = rhs.size;
}

// (Virtual) Destructor                                                                                                                                                       
template<typename T>
Operator<T>::~Operator() {}

// Assignment Operator                                                                                                                                                        
template<typename T>
Operator<T>& Operator<T>::operator=(const Operator<T>& rhs) {
	if (&rhs == this)
		return *this;

	this->elements = rhs.elements;
	this->size = rhs.size;

	return *this;
}
//TODO: once matched basis found break if for sure no redundancy exists    
// Addition of two Operators   
template<typename T>
Operator<T> Operator<T>::operator+(const Operator<T>& rhs) {
	std::deque<Projector<T>> new_elements;

	std::deque<Projector<T>> current_elements = this->elements;
	unsigned cur_elm_size = current_elements.size();

	std::deque<Projector<T>> rhs_elements = rhs.elements;
	
	for (unsigned i = cur_elm_size; i > 0; i--) {
		Projector<T> proj = current_elements.front();
		unsigned counter = 0;
		unsigned rhs_elm_size = rhs_elements.size();

		for (unsigned j = rhs_elm_size; j > 0; j--) {
			Projector<T> rhs_proj = rhs_elements.front();

			if (proj.match(rhs_proj))
			{
				counter += 1;
				proj.weight += rhs_proj.weight;
				//std::complex<double> complex_zero;

				//if (proj.weight != complex_zero)
				//{
					new_elements.push_back(proj);
					current_elements.pop_front();
					rhs_elements.pop_front();
					//break;
				//}
			}
			else
			{
				rhs_elements.push_back(rhs_proj);
				rhs_elements.pop_front();
			}
		}
		//iterated through all projectors in rhs and a match was not found
		if (counter == 0) {
			new_elements.push_back(proj);
			current_elements.pop_front();
		}
	}
	unsigned rhs_new_elm_size = rhs_elements.size();
	if (rhs_new_elm_size != 0)
	{
		for (unsigned j = rhs_new_elm_size; j > 0; j--) {
			new_elements.push_back(rhs_elements.front());
			rhs_elements.pop_front();
		}
	}

	this->size = rhs.size;

	Operator<T> result(this->size, new_elements);

	return result;
}

// Subtraction of this Operator and another  
//remove first element of current operator, remove first element of rhs deque, if they match, subtract weight and push_back to the result,
//if not matched push back the first elemnt of rhs to the end of deque and repeate the procedure, if iterataed through all elements of rhs and could not find a match, just
// add the first element of the current one to the result. repeat untill iterated through all elements of current deque (current deque be empty now), if rhs deque is not empty,
//all all remaining elements of rhs to the result
template<typename T>
Operator<T> Operator<T>::operator-(const Operator<T>& rhs) {
	std::deque<Projector<T>> new_elements;

	std::deque<Projector<T>> current_elements = this->elements;
	unsigned cur_elm_size = current_elements.size();

	std::deque<Projector<T>> rhs_elements = rhs.elements;

	for (unsigned i = cur_elm_size; i > 0; i--) {
		Projector<T> proj = current_elements.front();
		unsigned counter = 0;
		unsigned rhs_elm_size = rhs_elements.size();

		for (unsigned j = rhs_elm_size; j > 0; j--) {
			Projector<T> rhs_proj = rhs_elements.front();

			if (proj.match(rhs_proj))
			{
				counter += 1;
				proj.weight -= rhs_proj.weight;
				std::complex<double> complex_zero;

				//if (proj.weight != complex_zero)
				//{
					new_elements.push_back(proj);
					current_elements.pop_front();
					rhs_elements.pop_front();
					//break;
				//}
			}
			else
			{
				rhs_elements.push_back(rhs_proj);
				rhs_elements.pop_front();
			}
		}
		//iterated through all projectors in rhs and a match was not found
		if (counter == 0) {
			new_elements.push_back(proj);
			current_elements.pop_front();
		}
	}
	unsigned rhs_new_elm_size = rhs_elements.size();
	if (rhs_new_elm_size != 0)
	{
		for (unsigned j = rhs_new_elm_size; j > 0; j--) {
			new_elements.push_back(rhs_elements.front());
			rhs_elements.pop_front();
		}
	}

	this->size = rhs.size;

	Operator<T> result(this->size, new_elements);

	return result;
}

//TODO: use push and pop for deque
//// tensor product of this Operator and another                                                                                                                              
//template<typename T>
//Operator<T> Operator<T>::operator*(const Operator<T>& rhs) {
//	std::deque<Projector<T>> new_elements;
//	unsigned new_size = rhs.size * this->size;
//	
//	for (unsigned j = 0; j < rhs.get_elements(); j++) {
//		for (unsigned i = 0; i < this->get_elements(); i++) {
//			Projector<T> r = this->elements[i] * rhs.elements[j];
//			std::complex<double> compex_zero(0.0, 0.0);
//			if (r.weight != compex_zero)
//				new_elements.push_back(r);
//		}
//	}
//	Operator<T> result(new_size, new_elements);
//
//	return result;
//}

// tensor product of this Operator and another                                                                                                                              
template<typename T>
Operator<T> Operator<T>::operator*(const Operator<T>& rhs) {
	std::deque<Projector<T>> new_elements;
	unsigned new_size = rhs.size * this->size;

	for (unsigned j = 0; j < rhs.get_elements(); j++) {
		for (unsigned i = 0; i < this->get_elements(); i++) {
			Projector<T> r = this->elements[i] * rhs.elements[j];
			std::complex<double> complex_zero;
			if (r.weight != complex_zero)
				new_elements.push_back(r);
		}
	}

	Operator<T> result(new_size, new_elements);

	return result;
}

//Inner product of this Operator with another Operator
template<typename T>
Operator<T> Operator<T>::inner(const Operator<T>& rhs) {
	std::deque<Projector<T>> new_elements;

	for (unsigned j = 0; j < rhs.get_elements(); j++) 
	{
		for (unsigned i = 0; i < this->get_elements(); i++) 
		{
			Projector<T> r = this->elements[i].inner(rhs.elements[j]);

			std::complex<double> complex_zero;
			if (r.weight != complex_zero)
				new_elements.push_back(r);
		}
	}

	Operator<T> result(this->size, new_elements);

	return result;
}

// Return the transpose of this Operator                                                                                                                                       
template<typename T>
const Operator<T> Operator<T>::transpose() {
	std::deque<Projector<T>> new_elements;

	for (unsigned i = 0; i < Operator<T>::get_elements(); i++) {
		Projector<T> transposed = this->elements[i].transpose();
		//std::complex<double> compex_zero(0.0, 0.0);
		//if (r.weight != compex_zero)
		new_elements.push_back(transposed);
	}

	Operator<T> result(this->size, new_elements);

	return result;
}

// Return the transpose of this Operator                                                                                                                                       
template<typename T>
const Operator<T> Operator<T>::conjugate() {
	std::deque<Projector<T>> new_elements;

	for (unsigned i = 0; i < Operator<T>::get_elements(); i++) {
		Projector<T> conjugated = this->elements[i].conjugate();

		//if (conjugated.weight != 0)
		new_elements.push_back(conjugated);
	}

	Operator<T> result(this->size, new_elements);

	return result;
}

// Return the dagger of this Operator                                                                                                                                       
template<typename T>
const Operator<T> Operator<T>::dagger() {
	Operator<T> transposed = this->transpose();
	Operator<T> result = transposed.conjugate();

	return result;
}

// Operator/scalar multiplication                                                                                                                                               
template<typename T>
Operator<T> Operator<T>::operator*(const std::complex<double> rhs) {

	for (unsigned i = 0; i < this->get_elements(); i++) {
		this->elements[i].weight *= rhs;
	}

	return *this;
}

// Operator/scalar division                                                                                                                                                     
template<typename T>
Operator<T> Operator<T>::operator/(const std::complex<double> rhs) {

	for (unsigned i = 0; i < this->get_elements(); i++) {
		this->elements[i].weight /= rhs;
	}

	return *this;
}

// Access the individual elements                                                                                                                                             
template<typename T>
T& Operator<T>::operator()(const unsigned& index) {
	return this->elements[index];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& Operator<T>::operator()(const unsigned& index) const {
	return this->elements[index];
}

// Get the number of elements of Operator                                                                                                                                       
template<typename T>
unsigned Operator<T>::get_elements() const {
	return this->elements.size();
}

// Apply an operator to a state, a state is a combination of kets                                                                                                                                          
template<typename T>
State<T> Operator<T>::applyTo(const State<T>& rhs) {
	std::deque<Ket<T>> new_elements;
	// for each projector in this operator, projector.applyTo(ket)
	for (unsigned i = 0; i < this->get_elements(); i++) {
		Projector<T> proj = this->elements[i];
		std::deque<Ket<T>>::const_iterator it = rhs.elements.begin();
		while (it != rhs.elements.end()) {
			//for test
			Ket<T> rr = *it;
			Ket<T> r = proj.applyTo(*it);

			std::complex<double> zero;
			if (r.weight != zero) {
				new_elements.push_back(r);
			}
			++it;
		}
	}
	//TODO: new_elements may have equal ket_bases, to remove redundant kets and get final state, we need to add up all of the kets in new_elements
	
	std::deque<Ket<T>> result_kets;
	result_kets.push_back(new_elements.front());
	new_elements.pop_front();

	State<T> result(rhs.size, result_kets);

	while (new_elements.size() != 0) {
		result = result + new_elements.front();
		new_elements.pop_front();
	}

	return result;
}

template<typename T>
void Operator<T>::asMatrix(unsigned num_qubits, unsigned base, std::complex<double> *R) {
	//std::complex<double> zero;
	unsigned dim = this->size;
	unsigned N = pow(base, num_qubits);

	Utils::initZero<T>(N*N, R);
	//for (unsigned i = 0; i < (dim*dim); i++) {
	//	R[i] = zero;
	//}

	std::deque<Projector<T>>::const_iterator it = this->elements.begin();
	while (it != this->elements.end()) {
		Projector<T> p = *it;
		unsigned row = p.ket.asDecimal(base).first;
		unsigned col = p.bra.asDecimal(base).first;
		//printf("row: %u , col: %u\n", row, col);
		unsigned Rindex = col + row * N;
		R[col + row * N] = p.weight;
 		++it;
	}

	return;
}

template<typename T>
double * Operator<T>::asRealMatrix(unsigned base) {
	//std::complex<double> zero;
	//std::vector<std::complex<double>> result(this->size, zero);
	unsigned dim = this->size;
	double *R = (double *)mkl_malloc(dim*dim * sizeof(double), 64);

	for (unsigned i = 0; i < (dim*dim); i++) {
		R[i] = 0;
	}

	std::deque<Projector<T>>::const_iterator it = this->elements.begin();
	while (it != this->elements.end()) {
		Projector<T> p = *it;
		unsigned row = p.ket.asDecimal(base).first;
		unsigned col = p.bra.asDecimal(base).first;

		R[col + row * dim] = p.weight.real();
		++it;
	}

	return R;
}

template<typename T>
void Operator<T>::print() {		
	printf("Start printing of this Operator \n\n");
	for (unsigned i = 0; i < this->get_elements(); i++) {
		//TODO: remove variable after test
		Projector<T> proj = this->elements[i];
		proj.print();		
	}
	printf("End printing of this Operator \n\n");
}

#pragma endregion

#endif