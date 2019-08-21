#pragma once
// qubit is an abstract class, a qubit needs to know the neighbor qubits, it's label, and it's state

#ifndef __QUBIT_H
#define __QUBIT_H

#include <utility>
#include<bitset>
#include "dirac.h"
#include "utils.h"

using namespace std;
using namespace Utils;

#define PI 3.14159265

template <typename T> class Operator;
template <typename T> class State;

template <typename T> class Qubit {
public:
	
	Qtype type;
	pair<unsigned, unsigned> index;

	//vector<pair<unsigned, unsigned>> entanglment_list; // initially this is an empty vector, during the evolution the elements will be added or removed. 
	// Maybe The network needs to know entangled indexes and Qubit only need to know if it is entangled!
	bool isEntangled;

	// Qubit needs to know its state and network should be able to write the state of qubits, also to read the state of qubits at any time
	// if entangled the network handles the state of entangled qubits altogether.
	// When constructing the product state, only qubits with isEntangled = false participate in tensor product, then Network will tensor the result with entangled state

	//state of a qubit is in the form aKet(0)+bKet(1)+cKet(2) + ...+ zKet(r) r is the number of atom levels, it can be represented by a vector of probability aplitudes
	// probAmplitudes[r] equals the weight of the r'th basis ket, so the vector always has r+1 elements (r is max atom level considered)
	//vector<double> probAmplitudes; //The summation of pow(probAmplitudes[i], 2) of all elements (i= 0, 1, ..., r) equals 1 
	State<T> state; //TODO: initialize to ground state

	Qubit();
	Qubit(Qtype q_type, pair<unsigned, unsigned> _index, bool isEntangled = false);
	virtual ~Qubit();

	//In the beginning, and then after each gate operation, the final state is updated using this method
	//while the state transtion history is saved to a file for plotting purpose
	void setState(State<T> state);

	State<T> getState();

	virtual Operator<T> get_hamiltonian() = 0;

	virtual Operator<T> get_effective_hamiltonian() = 0;

	virtual Operator<T> get_hamiltonian_drive(double drive_amplitude, double Wd) = 0;

	virtual Operator<T> get_effective_hamiltonian_drive(double drive_amplitude, double Wd) = 0;
};

#include "qubit.cpp"

#endif