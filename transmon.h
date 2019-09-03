#pragma once

#ifndef __TRANSMON_H
#define __TRANSMON_H

#include <math.h>
#include <utility>
#include<bitset>
#include "dirac.h"
#include "utils.h"

using namespace std;
using namespace Utils;

#define PI 3.14159265
// A transmon inherits from a qubit, transmon knows its characteristics such as transition freq, cavity coupling, etc, and can implement Hamiltonians
//

template <typename T> class Qubit;
template <typename T> class Operator;

template <typename T> class Transmon : public Qubit<T>
{
public:
	double e = pow(1.60217662, -19); //electron charge	
	double h = pow(6.62607004, -34); //planck's constant (m^2kg/s) or (J.s) we consider 1 for easier simulation
	double phi0 = h / (2 * e); //superconducting flux quantum

	double Ej1 = 2 * PI * 4500; //MJ
	double Ej2 = 2 * PI * 4500;
	double Ej_sigma = Ej1 + Ej2;
	double Cj = 0.009 * pow(10, -6); //Josephson capacitance in microFarad

	double Cg = 0.08 * pow(10, -6); //additional capacitance in microFarad
	
	//double Vg; //gate voltage
	//double ng = Cg * Vg / (2 * e); //the charge offset on the island can be controlled by Vg, in sweet-spot of coherence ng=1/2

	double C = Cg + Cj;
	double beta = Cg / C;
	//TOD: define this
	double V0rms;

	//int n; //number of cooper pairs that crossed the junction
	//double q = -n * 2 * e;//discrete charge

	double Ec = pow(e, 2) / (2 * C); //charging energy
	double phi;//flux through the josephson junction is a control parameter to change the transition frequency

	double d = 0.1; // d = (Ej2 - Ej1) / Ej_sigma;  junction asymmetry

	double Ej = Ej_sigma * cos(PI * phi / phi0) * sqrt(1 + pow(d, 2) * pow(tan(PI * phi / phi0), 2)); //josephson energy for SQUID 

	double g = pow((Ej / 2 * Ec), 1.0 / 4) * e * beta * V0rms;//atom-cavity coupling

	double w0 = sqrt(8 * Ec * Ej);
	double delta = -Ec; //anharmonicity
	double w = w0 + delta;

	unsigned r = 2; //max considered atom levels
	//Qtype type;
	//pair<unsigned, unsigned> index;

	//vector<pair<unsigned, unsigned>> entanglment_list; // initially this is an empty vector, during the evolution the elements will be added or removed. 
	// Maybe The network needs to know entangled indexes and Qubit only need to know if it is entangled!
	//bool isEntangled;

	// Qubit needs to know its state and network should be able to write the state of qubits, also to read the state of qubits at any time
	// if entangled the network handles the state of entangled qubits altogether.
	// When constructing the product state, only qubits with isEntangled = false participate in tensor product, then Network will tensor the result with entangled state

	//state of a qubit is in the form aKet(0)+bKet(1)+cKet(2) + ...+ zKet(r) r is the number of atom levels, it can be represented by a vector of probability aplitudes
	// probAmplitudes[r] equals the weight of the r'th basis ket, so the vector always has r+1 elements (r is max atom level considered)
	//vector<double> probAmplitudes; //The summation of pow(probAmplitudes[i], 2) of all elements (i= 0, 1, ..., r) equals 1 

	T photon_number; // average photon number
	T atom_level; // max atom level
	double Wr; // Cavity_frequency

	Transmon();
	Transmon(Qtype q_type, pair<unsigned, unsigned> _index, T& _photon_number, T& _atom_level, double _cavity_frequency, double transition_freq, double qubit_cavity_coupling, double anharmonicity, bool isEntangled = false);// : Qubit<T>(q_type, _index, isEntangled);
	virtual ~Transmon();

	//In the beginning, and then after each gate operation, the final state is updated using this method
	//while the state transtion history is saved to a file for plotting purpose
	//void updateState(State<T> state);

	//cosidered atom level as unsigned int
	double get_trans_freq(const T& j, double w);

	Operator<T> get_trans_freq_operator(const T& max_atom_level, double w0);

	double get_mu(const T& j, double Wr, double w);

	double get_J(const T& j1, const T& j2, double Wr, double W1, double W2, double delta1, double delta2, double g1, double g2);

	Operator<T> Hcavity(double _Wr);

	Operator<T> Hbare();

	Operator<T> Hcoupling();


	Operator<T> get_hamiltonian();

	Operator<T> get_effective_hamiltonian();


	Operator<T> get_hamiltonian_drive(double drive_amplitude, double Wd);

	Operator<T> get_effective_hamiltonian_drive(double drive_amplitude, double Wd);

	std::pair<Operator<T>, double> get_Swap_gate(Transmon<T> target, double drive_amplitude, double Wd);

	std::pair<Operator<T>, double> Transmon<T>::get_X_Y_gate(double drive_x_amplitude, double drive_y_amplitude, double Wd, double x_theta, double y_theta);
	
	Operator<T> Transmon<T>::get_X_or_Y_gate(double drive_amplitude, double Wd, double x_theta, double y_theta);

	Operator<T> Transmon<T>::get_Cphase_gate_freq_pulse(Transmon<T> target, double w_l, double w_r);
};

#include "transmon.cpp"

#endif
