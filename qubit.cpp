#ifndef __QUBIT_CPP
#define __QUBIT_CPP

#include "qubit.h"

// Default Constructor    
template<typename T>
Qubit<T>::Qubit() {}

template<typename T>
Qubit<T>::Qubit(Qtype q_type, pair<unsigned, unsigned> _index, bool _isEntangled = false){
	type = q_type;
	index = _index;
	isEntangled = _isEntangled;
}

// (Virtual) Destructor 
template<typename T>
Qubit<T>::~Qubit() {}

template<typename T>
void Qubit<T>::setState(State<T> _state) {
	this->state = _state;

	return;
}

template<typename T>
State<T> Qubit<T>::getState() {
	State<T> state = this->state;

	return state;
}

#endif
