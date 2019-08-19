//https://baptiste-wicht.com/posts/2012/12/cpp-benchmark-vector-list-deque.html
#ifndef __DIRAC_CPP
#define __DIRAC_CPP

#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>
#include <string>
#include <stdexcept>

#include "dirac.h"

//template <typename T> class State;

//template <typename T> class Projector {
//public:
//	Bra<T> bra;
//	Ket<T> ket;
//	std::complex<double> weight;
//
//	// Default Constructor                                                                                                                                                      
//	Projector();
//
//	Projector(Bra<T> bra, Ket<T> ket, std::complex<double> weight);
//
//	Projector(std::vector<T> bra_elements, std::vector<T> ket_elements, std::complex<double> weight);
//
//	// Copy Constructor    
//	Projector(const Projector<T>& rhs);
//	virtual ~Projector();
//
//	// Operator overloading, for "standard" mathematical dirac operations                                                                                                                                                          
//	Projector<T>& operator=(const Projector<T>& rhs);
//
//	// Dirac mathematical operations                                                                                                                                                                                               
//	std::deque<Projector<T>> operator+(const Projector<T>& rhs);
//	std::deque<Projector<T>> operator-(const Projector<T>& rhs);
//	Projector<T> operator*(const Projector<T>& rhs);
//
//	Projector<T> inner(const Projector<T>& rhs);
//	//Projector<T> innerleft(const Projector<T>& lhs);
//	Projector<T> transpose();
//	Projector<T> conjugate();
//
//	// Bra/scalar multiplication
//	Projector<T> operator*(const T& rhs);
//	// Bra/scalar division                                                                                                                                                     
//	Projector<T> operator/(const T& rhs);
//
//	Ket<T> applyTo(const Ket<T>& rhs);
//
//	// Ket/vector operations                                                                                                                                                                                                     
//	//std::vector<T> operator*(const std::vector<T>& rhs);
//	//std::vector<T> diag_vec();
//	bool matchKet(const Projector<T> rhs);
//	bool matchBra(const Projector<T> rhs);
//	bool match(const Projector<T> rhs);
//};

//template <typename T> class Operator {
//public:
//	std::deque<Projector<T>> elements;
//	unsigned size; //If the Operator is like a n by n matrix, the size is n
//
//	Operator();
//	Operator(unsigned size, const std::deque<Projector<T>> elements);
//	Operator(const Operator<T>& rhs);
//	virtual ~Operator();
//
//	// Operator overloading, for "standard" mathematical dirac operations                                                                                                                                                          
//	Operator<T>& operator=(const Operator<T>& rhs);
//
//	// Dirac mathematical operations                                                                                                                                                                                               
//	Operator<T> operator+(const Operator<T>& rhs);
//	Operator<T> operator-(const Operator<T>& rhs);
//	Operator<T> operator*(const Operator<T>& rhs);
//
//	// Inner product of this Operator with a Operator
//	Operator<T> inner(const Operator<T>& rhs);
//	Operator<T> transpose();
//	Operator<T> conjugate();
//	Operator<T> dagger();
//
//	//TODO: Commutator
//
//	//Projector<T> outer(const Ket<T>& lhs);
//	//Projector<T&> outer(const std::vector<T&> ket_elements, unsigned ket_weight, T&);
//
//	// Bra/scalar multiplication
//	Operator<T> operator*(const std::complex<float> rhs);
//	// Bra/scalar division                                                                                                                                                     
//	Operator<T> operator/(const std::complex<float> rhs);
//
//	// Bra/vector operations                                                                                                                                                                                                     
//	//std::vector<T> operator*(const std::vector<T>& rhs);
//	//std::vector<T> diag_vec();
//
//	// Access the individual element                                                                                                                                                                                               
//	T& operator()(const unsigned& element);
//	const T& operator()(const unsigned& element) const;
//
//	// Get the number of elements of operator                                                                                                                                                                                           
//	unsigned get_elements() const;
//
//	State<T> applyTo(const State<T>& rhs);
//};

//template <typename T> class State {
//public:
//	std::deque<Ket<T>> elements;
//	unsigned size;
//
//	State();
//	State(unsigned size, const std::deque<Ket<T>> elements);
//	State(const State<T>& rhs);
//	virtual ~State();
//
//	// Operator overloading, for "standard" mathematical dirac operations                                                                                                                                                          
//	State<T>& operator=(const State<T>& rhs);
//
//	// State/scalar multiplication
//	State<T> operator*(const T& rhs);
//	// State/scalar division                                                                                                                                                     
//	State<T> operator/(const T& rhs);
//
//	// Access the individual element                                                                                                                                                                                               
//	std::deque<Ket<T>>& operator()(const unsigned& element);
//	const std::deque<Ket<T>>& operator()(const unsigned& element) const;
//
//	// Get the number of elements of operator                                                                                                                                                                                           
//	unsigned get_elements() const;
//};

//#pragma region Projector
//
//// Default Constructor                                                                                                                                                      
//template<typename T>
//Projector<T>::Projector() {}
//
//// Parameter Constructor 
//template<typename T>
//Projector<T>::Projector(Bra<T> _bra, Ket<T> _ket, std::complex<double> _weight) {
//	bra = _bra;
//	ket = _ket;
//	weight = _weight;
//}
//
//template<typename T>
//Projector<T>::Projector(std::vector<T> _bra_elements, std::vector<T> _ket_elements, std::complex<double> _weight) {
//	Bra<T> _bra(1, _bra_elements, 1.0);
//	Ket<T> _ket(1, _ket_elements, 1.0);
//	bra = _bra;
//	ket = _ket;
//	weight = _weight;
//}
//
//// Copy Constructor                                                                                                                                                           
//template<typename T>
//Projector<T>::Projector(const Projector<T>& rhs) {
//	bra = rhs.bra;
//	ket = rhs.ket;
//	weight = rhs.weight;
//}
//
//
//// (Virtual) Destructor                                                                                                                                                       
//template<typename T>
//Projector<T>::~Projector() {}
//
//// Assignment Operator                                                                                                                                                        
//template<typename T>
//Projector<T>& Projector<T>::operator=(const Projector<T>& rhs) {
//	if (&rhs == this)
//		return *this;
//
//	this->bra = rhs.bra;
//	this->ket = rhs.ket;
//	this->weight = rhs.weight;
//
//	return *this;
//}
//
//// Addition of two Projectors                                                                                                                                                 
//template<typename T>
//std::deque<Projector<T>> Projector<T>::operator+(const Projector<T>& rhs) {
//	std::deque<Projector<T>> result;
//
//	if (this->match(rhs)) {
//		this->weight += rhs.weight;		
//	}
//	else
//	{
//		result.push_back(rhs);
//	}
//
//	result.push_back(this);
//
//	return result;
//}
//#pragma warning: TODO: update subtraction
//
//// Subtraction of this Projector and another                                                                                                                                     
//template<typename T>
//std::deque<Projector<T>> Projector<T>::operator-(const Projector<T>& rhs) {
//	std::deque<Projector<T>> result;
//	if (this->match(rhs)) {
//		this->weight -= rhs.weight;
//	}
//	else
//	{
//		result.push_back(rhs);
//	}
//
//	result.push_back(this);
//
//	return result;
//}
//
//// Left tensor product of this Projector and another                                                                                                                              
//template<typename T>
//Projector<T> Projector<T>::operator*(const Projector<T>& rhs) {
//	Bra<T> result_bra = this->bra * rhs.bra;
//	Ket<T> result_ket = this->ket * rhs.ket;
//	std::complex<double> result_weight = this->weight * rhs.weight;
//
//	Projector<T> result(result_bra, result_ket, result_weight);
//
//	return result;
//}
//
//// Inner product of this Projector with another Projector applied to right	
//// If bra_ket_result is zero, the projector weight will be zero.
//template<typename T>
//Projector<T> Projector<T>::inner(const Projector<T>& rhs) {
//	std::complex<double> bra_ket_result = this->bra.inner(rhs.ket);
//	std::complex<double> result_weight = bra_ket_result * this->weight * rhs.weight;
//
//	Projector<T> result(rhs.bra, this->ket, result_weight);
//
//	return result;
//}
//
//// Return the transpose of this projector                                                                                                                                       
//template<typename T>
//Projector<T> Projector<T>::transpose() {
//	Projector<T> result(this->ket.elements, this->bra.elements, this->weight);
//
//	return result;
//}
//
//// Return the conjugate of this projector                                                                                                                                       
//template<typename T>
//Projector<T> Projector<T>::conjugate() {
//	this->weight = std::conj(this->weight);
//
//	return *this;
//}
//
////TODO: implement dagger
//
//// Projector/scalar multiplication                                                                                                                                               
//template<typename T>
//Projector<T> Projector<T>::operator*(const T& rhs) {
//	this->weight *= rhs
//
//		return *this;
//}
//
//// Projector/scalar division                                                                                                                                                     
//template<typename T>
//Projector<T> Projector<T>::operator/(const T& rhs) {
//	this->weight /= rhs
//
//		return *this;
//}
//
//// Apply a Projector to a ket                                                                                                                                            
//template<typename T>
//Ket<T> Projector<T>::applyTo(const Ket<T>& rhs) {
//	double result_weight = this->bra.inner(rhs) * this->weight * rhs.weight;
//
//	Ket<T> result(rhs.length, this->ket, result_weight);
//
//	return result;
//}
//
//// Returns true if this projecr ket basis matches the rhs's ket basis                                                                                                                                          
//template<typename T>
//bool Projector<T>::matchKet(const Projector<T> rhs) {	
//
//	return this->ket.elements == rhs.ket.elements;
//}
//
//// Returns true if this projecr bra basis matches the rhs's bra basis                                                                                                                                          
//template<typename T>
//bool Projector<T>::matchBra(const Projector<T> rhs) {
//
//	return this->bra.elements == rhs.bra.elements;
//}
//
//// Returns true if this projecr bra basis matches the rhs's bra basis                                                                                                                                          
//template<typename T>
//bool Projector<T>::match(const Projector<T> rhs) {
//
//	return matchKet(rhs) && matchBra(rhs);
//}
//
//#pragma endregion

// Change data structure of Operaqtor to map of key(projector.bra), value (projector)
//#pragma region Operator
//
//// Default Constructor                                                                                                                                                      
//template<typename T>
//Operator<T>::Operator() {}
//
//// Parameter Constructor                                                                                                                                                      
//template<typename T>
//Operator<T>::Operator(unsigned _size, const std::deque<Projector<T>> _elements) {
//	elements.resize(_elements.size());
//	elements = _elements;
//	size = _size;
//}
//
//// Copy Constructor                                                                                                                                                           
//template<typename T>
//Operator<T>::Operator(const Operator<T>& rhs) {
//	elements = rhs.elements;
//	size = rhs.size;
//}
//
//// (Virtual) Destructor                                                                                                                                                       
//template<typename T>
//Operator<T>::~Operator() {}
//
//// Assignment Operator                                                                                                                                                        
//template<typename T>
//Operator<T>& Operator<T>::operator=(const Operator<T>& rhs) {
//	if (&rhs == this)
//		return *this;
//
//	this->elements = rhs.elements;
//	this->size = rhs.size;
//
//	return *this;
//}
////TODO: once matched basis found break if for sure no redundancy exists    
//// Addition of two Operators   
//template<typename T>
//Operator<T> Operator<T>::operator+(const Operator<T>& rhs) {
//	std::deque<Projector<T>> new_elements;
//
//	std::deque<Projector<T>> current_elements = this->elements;
//	unsigned cur_size = current_elements.size();
//
//	std::deque<Projector<T>> rhs_elements = rhs.elements;
//	unsigned rhs_size = rhs_elements.size();
//
//	for (unsigned i = size; i > 0; i--) {
//		Projector<T> proj = current_elements.front();
//		unsigned counter = 0;
//
//		for (unsigned j = rhs_size; j > 0; j--) {
//			Projector<T> rhs_proj = rhs_elements.front();
//
//			if (proj.match(rhs_proj))
//			{
//				counter += 1;
//				proj.weight += rhs_proj.weight;
//
//				if (proj.weight != 0.0)
//				{
//					new_elements.push_back(proj);
//					current_elements.pop_front();
//					rhs_elements.pop_front();
//				}
//			}
//			else
//			{
//				rhs_elements.push_back(rhs_proj);
//				rhs_elements.pop_front();
//			}
//		}
//		//iterated through all projectors in rhs and a match was not found
//		if (counter == 0) {
//			new_elements.push_back(proj);
//			current_elements.pop_front();
//		}
//	}
//	unsigned rhs_new_size = rhs_elements.size();
//	if (rhs_new_size != 0)
//	{
//		for (unsigned j = rhs_new_size; j > 0; j--) {
//			new_elements.push_back(rhs_elements.front());
//			rhs_elements.pop_front();
//		}
//	}
//	Operator<T> result(new_elements.size(), new_elements);
//
//	return result;
//}
//
//// Subtraction of this Operator and another  
////remove first element of current operator, remove first element of rhs deque, if they match, subtract weight and push_back to the result,
////if not matched push back the first elemnt of rhs to the end of deque and repeate the procedure, if iterataed through all elements of rhs and could not find a match, just
//// add the first element of the current one to the result. repeat untill iterated through all elements of current deque (current deque be empty now), if rhs deque is not empty,
////all all remaining elements of rhs to the result
//template<typename T>
//Operator<T> Operator<T>::operator-(const Operator<T>& rhs) {
//	std::deque<Projector<T>> new_elements;
//
//	std::deque<Projector<T>> current_elements = this->elements;
//	unsigned cur_size = current_elements.size();
//
//	std::deque<Projector<T>> rhs_elements = rhs.elements;
//	unsigned rhs_size = rhs_elements.size();
//
//	for (unsigned i = size; i > 0; i--) {
//		Projector<T> proj = current_elements.front();
//		unsigned counter = 0;
//
//		for (unsigned j = rhs_size; j > 0; j--) {
//			Projector<T> rhs_proj = rhs_elements.front();
//
//			if (proj.match(rhs_proj))
//			{
//				counter += 1;
//				proj.weight -= rhs_proj.weight;
//
//				if (proj.weight != 0.0)
//				{
//					new_elements.push_back(proj);
//					current_elements.pop_front();
//					rhs_elements.pop_front();
//				}
//			}
//			else
//			{
//				rhs_elements.push_back(rhs_proj);
//				rhs_elements.pop_front();
//			}
//		}
//		//iterated through all projectors in rhs and a match was not found
//		if (counter == 0) {
//			new_elements.push_back(proj);
//			current_elements.pop_front();
//		}
//	}
//	unsigned rhs_new_size = rhs_elements.size();
//	if (rhs_new_size != 0)
//	{
//		for (unsigned j = rhs_new_size; j > 0; j--) {
//			new_elements.push_back(rhs_elements.front());
//			rhs_elements.pop_front();
//		}
//	}
//	Operator<T> result(new_elements.size(), new_elements);
//
//	return result;
//}
//
//// tensor product of this Operator and another                                                                                                                              
//template<typename T>
//Operator<T> Operator<T>::operator*(const Operator<T>& rhs) {
//	std::deque<Projector<T>> new_elements;
//	unsigned new_size = rhs.size * this->size;
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
//
////Inner product of this Operator with another Operator
//template<typename T>
//Operator<T> Operator<T>::inner(const Operator<T>& rhs) {
//	std::deque<Projector<T>> new_elements;
//	for (unsigned j = 0; j < rhs.get_elements(); j++) {
//		for (unsigned i = 0; i < this->get_elements(); i++) {
//			Projector<T> r = this->elements[i].inner(rhs.elements[j]);
//
//			std::complex<double> compex_zero(0.0, 0.0);
//			if (r.weight != compex_zero)
//				new_elements.push_back(r);
//		}
//	}
//	Operator<T> result(this->size, new_elements);
//
//	return result;
//}
//
//// Return the transpose of this Operator                                                                                                                                       
//template<typename T>
//Operator<T> Operator<T>::transpose() {
//	std::deque<Projector<T>> new_elements;
//	for (unsigned i = 0; i < Operator<T>::get_elements(); i++) {
//		Projector<T> transposed = this->elements[i].transpose();
//		//std::complex<double> compex_zero(0.0, 0.0);
//		//if (r.weight != compex_zero)
//		new_elements.push_back(transposed);
//	}
//	Operator<T> result(this->size, new_elements);
//
//	return result;
//}
//
//// Return the transpose of this Operator                                                                                                                                       
//template<typename T>
//Operator<T> Operator<T>::conjugate() {
//	std::deque<Projector<T>> new_elements;
//	for (unsigned i = 0; i < Operator<T>::get_elements(); i++) {
//		Projector<T> conjugated = this->elements[i].conjugate();
//
//		//if (conjugated.weight != 0)
//		new_elements.push_back(conjugated);
//	}
//	Operator<T> result(this->size, new_elements);
//
//	return result;
//}
//
//// Return the dagger of this Operator                                                                                                                                       
//template<typename T>
//Operator<T> Operator<T>::dagger() {
//	Operator<T> transposed =this->transpose();
//	Operator<T> result = transposed.conjugate();
//
//	return result;
//}
//
//// Operator/scalar multiplication                                                                                                                                               
//template<typename T>
//Operator<T> Operator<T>::operator*(const std::complex<float> rhs) {
//	for (unsigned i = 0; i < this->get_elements(); i++) {
//		this->elements[i].weight *= rhs;
//	}
//
//	return *this;
//}
//
//// Operator/scalar division                                                                                                                                                     
//template<typename T>
//Operator<T> Operator<T>::operator/(const std::complex<float> rhs) {
//	for (unsigned i = 0; i < this->get_elements(); i++) {
//		this->elements[i].weight /= rhs;
//	}
//
//	return *this;
//}
//
//// Access the individual elements                                                                                                                                             
//template<typename T>
//T& Operator<T>::operator()(const unsigned& index) {
//	return this->elements[index];
//}
//
//// Access the individual elements (const)                                                                                                                                     
//template<typename T>
//const T& Operator<T>::operator()(const unsigned& index) const {
//	return this->elements[index];
//}
//
//// Get the number of elements of Operator                                                                                                                                       
//template<typename T>
//unsigned Operator<T>::get_elements() const {
//	return this->elements.size();
//}
//
//// Apply an operator to a state                                                                                                                                            
//template<typename T>
//State<T> Operator<T>::applyTo(const State<T>& rhs) {
//	std::deque<Ket<T>> new_elements;
//	// for each projector in this operator, projector.applyTo(ket)
//	for (unsigned i = 0; i < this->get_elements(); i++) {
//		Projector<T> proj = this->elements[i];
//		
//		for (std::deque<Ket<T>>::iterator it = rhs.elements.begin(); it != rhs.elements.end(); ++it) {
//			Ket<T> r = proj.applyTo(*it);
//			
//			if (r.weight != 0) {
//				new_elements.push_back(r);
//			}
//		}
//	}
//	//TODO: new_elements may have equal ket_bases, to remove redundant kets and get final state, we need to add up all of the kets in new_elements
//	//Ket<T> result = new_elements.front();
//	//new_elements.pop_front();
//
//	//while (new_elements.size() != 0) {
//	//	result = result + new_elements.front();
//	//	new_elements.pop_front();
//	//}
//
//	State<T> result_state(new_elements.size, new_elements);
//
//	return result_state;
//}
//
//#pragma endregion

//#pragma region State
//
//// Default Constructor                                                                                                                                                      
//template<typename T>
//State<T>::State() {}
//
//// Parameter Constructor                                                                                                                                                      
//template<typename T>
//State<T>::State(unsigned _size, const std::deque<Ket<T>> _elements) {
//	elements.resize(_elements.size());
//	elements = _elements;
//	size = _size;
//}
//
//// Copy Constructor                                                                                                                                                           
//template<typename T>
//State<T>::State(const State<T>& rhs) {
//	elements = rhs.elements;
//}
//
//// (Virtual) Destructor                                                                                                                                                       
//template<typename T>
//State<T>::~State() {}
//
//// Assignment Operator                                                                                                                                                        
//template<typename T>
//State<T>& State<T>::operator=(const State<T>& rhs) {
//	if (&rhs == this)
//		return *this;
//
//	this->elements = rhs.elements;
//	this->size = rhs.size;
//	return *this;
//}
//
//// Operator/scalar multiplication                                                                                                                                               
//template<typename T>
//State<T> State<T>::operator*(const T& rhs) {
//	for (unsigned i = 0; i < Operator<T>::get_elements(); i++) {
//		this->elements(i).weight *= rhs
//	}
//
//	return *this;
//}
//
//// Operator/scalar division                                                                                                                                                     
//template<typename T>
//State<T> State<T>::operator/(const T& rhs) {
//	for (unsigned i = 0; i < Operator<T>::get_elements(); i++) {
//		this->elements(i).weight /= rhs
//	}
//
//	return *this;
//}
//
//// Access the individual elements                                                                                                                                             
//template<typename T>
//std::deque<Ket<T>>& State<T>::operator()(const unsigned& index) {
//
//	return this->elements[index];
//}
//
//// Access the individual elements (const)                                                                                                                                     
//template<typename T>
//const std::deque<Ket<T>>& State<T>::operator()(const unsigned& index) const {
//
//	return this->elements[index];
//}
//
//// Get the number of elements of Operator                                                                                                                                       
//template<typename T>
//unsigned State<T>::get_elements() const {
//
//	return this->elements.size();
//}
//
//#pragma endregion

#endif