#pragma once
// Network is an array of qubits, a network can update the state of its qubit, or perform gate operations on them, for two qubits gate check if the qubits are neighbor
// Network can mark qubits with types such as Data, X stabilizer, Z stabilizer, and O which means it is considered don't care for the purpose of calculation
#pragma once

#ifndef __NETWORK_H
#define __NETWORK_H

#include "dirac.h"

#include  <vector>
#include <complex>
#include <map>
#include <iterator>
#include <deque>

template <typename T> class Qubit;

template <typename T> class Network {
public:
	std::vector<Qubit<T>> qubits;

	Network();
	Network(const std::vector<Qubit<T>> qubits);
	virtual ~Network();

};

#include "network.cpp"

#endif