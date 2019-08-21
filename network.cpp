#ifndef __NETWORK_CPP
#define __NETWORK_CPP

#include "network.h"

#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

#pragma region Network

// Default Constructor                                                                                                                                                      
template<typename T>
Network<T>::Network() {}

// Parameter Constructor                                                                                                                                                      
template<typename T>
Network<T>::Network(const std::vector<Qubit<T>> _qubits){
	qubits = _qubits;
}

// (Virtual) Destructor                                                                                                                                                       
template<typename T>
Network<T>::~Network() {}

#pragma endregion

#endif