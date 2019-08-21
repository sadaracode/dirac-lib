#ifndef __TRANSMON_CPP
#define __TRANSMON_CPP

#include "transmon.h"

// Default Constructor    
template<typename T>
Transmon<T>::Transmon() {}

//add parameters as one wrapped class
template<typename T>
Transmon<T>::Transmon(Qtype q_type, pair<unsigned, unsigned> _index, T& _photon_number, T& _atom_level, double _cavity_frequency, double transition_freq, double qubit_cavity_coupling, double anharmonicity, bool _isEntangled = false)
			: Qubit<T>(q_type, _index, isEntangled) 
{
	photon_number = _photon_number;
	atom_level = _atom_level;
	Wr = _cavity_frequency;
	w = transition_freq;
	g = qubit_cavity_coupling;	
	delta = anharmonicity;
}

// (Virtual) Destructor 
template<typename T>
Transmon<T>::~Transmon() {}

template<typename T>
double Transmon<T>::get_trans_freq(const T& j, double w) {

	double result = ((w - (this->delta / 2)) * j) + ((this->delta / 2) * pow(j, 2));
	return result;
}

template<typename T>
Operator<T> Transmon<T>::get_trans_freq_operator(const T& max_atom_level, double w0) 
{
	Operator<T> photon_number = creation(max_atom_level).inner(annihilation(max_atom_level));
	Operator<T> H = (photon_number * (w0 - (this->delta / 2))) + (photon_number.inner(photon_number) * (this->delta) / 2);

	return H;
}

template<typename T>
double Transmon<T>::get_mu(const T& j, double Wr, double w) {
	double result = j / (Wr - w - (this->delta * (j - 1)));

	return result;
}

template<typename T>
double Transmon<T>::get_J(const T& j1, const T& j2, double Wr, double W1, double W2, double delta1, double delta2, double g1, double g2) {
	double result = g1 * g2 * (W1 + W2 + (delta1*j1) + (delta2*j2 - 2 * Wr)) / (2 * (Wr - W2 - delta2 * j2)*(Wr - W1 - delta1 * j1));

	return result;
}

template<typename T>
Operator<T> Transmon<T>::Hcavity(double _Wr) {
	Operator<T> a = annihilation(this->photon_number);
	Operator<T> a_dagger = a.dagger();
	T n = 1;
	//Operator<T> H =  Utils::pauliI<T>() * n * _Wr;
	Operator<T> H = a_dagger.inner(a) * _Wr;
	return H;
}

template<typename T>
Operator<T> Transmon<T>::Hbare() {
	std::deque<Projector<T>> level_projs;	
	T j = this->atom_level;

	for (T k = 0; k < j+1; k++) {
		double freq = this->get_trans_freq(k, this->w);
		if (freq != 0.0) 
		{
			Projector<T> proj({ k }, { k }, freq);
			level_projs.push_back(proj);
		}
	}
	Operator<T> H(j, level_projs);

	return H;
}

template<typename T>
Operator<T> Transmon<T>::Hcoupling() {
	Operator<T> a = annihilation(this->photon_number);
	Operator<T> a_dagger = a.dagger();

	Operator<T> c = annihilation(this->atom_level);
	Operator<T> c_dagger = c.dagger();

	Operator<T> coupling = a_dagger.inner(c) + a.inner(c_dagger);

	Operator<T> H = coupling * this->g;
	return H;
}

template<typename T>
Operator<T> Transmon<T>::get_hamiltonian() {
	Operator<T> Htemp = Hcavity(this->Wr) + Hbare();
	Operator<T> H = Htemp + Hcoupling();

	return H;
}

template<typename T>
Operator<T> Transmon<T>::get_effective_hamiltonian() {
	T n = this->photon_number;
	T j = this->atom_level;

	Operator<T> a = annihilation(n);
	Operator<T> a_dagger = a.dagger();
	a.print();
	a_dagger.print();

	std::deque<Projector<T>>kay_projs;
	//TODO: Check this
	//for (T k = 0; k < j ; k++)
	for (T k = 0; k < j+1; k++) 
	{
		//TODO: check this
		double kay = pow(g, 2)* (this->get_mu(k+1, this->Wr, this->w) - this->get_mu(k, this->Wr, this->w));
		if (kay != 0.0)
		{
			Projector<T> proj({ k }, { k }, kay);
			kay_projs.push_back(proj);
		}
	}
	
	Operator<T> Hkay(j+1, kay_projs);
	Hkay.print();

	std::deque<Projector<T>>trans_projs;

	for (T k = 0; k < j + 1; k++) {
		double Wbar = this->get_trans_freq(k, this->w) - (pow(g, 2)* this->get_mu(k, this->Wr, this->w));
		if (Wbar != 0.0)
		{
			Projector<T> proj({ k }, { k }, Wbar);
			trans_projs.push_back(proj);
		}
	}
	Operator<T> Htrans_eff(j+1, trans_projs);
	Htrans_eff.print();
	
	//Operator<T> H = this->Hcavity(this->Wr) + Hkay.inner(a_dagger.inner(a)) + Htrans_eff;
	//Operator<T> H = this->Hcavity(this->Wr) + Hkay.inner(Utils::pauliI<T>()) + Htrans_eff;
	Operator<T> H = this->Hcavity(this->Wr) + Utils::pauliI<T>().inner(Hkay) + Htrans_eff;
	H.print();

	return H;
}

//qubit approximation: two atom levels qubit, j=1, n=1
template<typename T>
Operator<T> Transmon<T>::get_hamiltonian_drive(double drive_amplitude, double Wd) {
	double delta_r = abs(this->Wr - Wd);
	double delta_q = abs(this->w - Wd);
	double _delta = abs(Wr - this->w);

	Operator<T> Hc = Hcavity(delta_r);

	Operator<T> Ht = pauliZ<T>() * (-delta_q / 2);

	std::complex<double> omega = -2 * (g / _delta) * drive_amplitude;
	Operator<T> Hd = (pauliMinus<T>() * std::conj(omega) + pauliPlus<T>() * omega) * (1 / 2);

	Operator<T> H = Hc + Ht + this->Hcoupling() + Hd;

	return H;
}

//qubit approximation: two atom levels qubit, j=1, n=1
template<typename T>
Operator<T> Transmon<T>::get_effective_hamiltonian_drive(double drive_amplitude, double Wd) {
	double delta_r = abs(this->Wr - Wd);
	double delta_q = abs(this->w - Wd);
	double _delta = abs(Wr - this->w);

	Operator<T> a = annihilation(this->photon_number);
	Operator<T> a_dagger = a.dagger();

	Operator<T> Hc = Hcavity(delta_r);

	Operator<T> Hz = pauliZ<T>().inner(a_dagger.inner(a)) * (pow(g, 2)/_delta) + pauliZ<T>() * (pow(g, 2) / (2* _delta));
	//Operator<T> Hz = pauliZ<T>().inner(Utils::pauliI<T>()) * ((pow(g, 2) / _delta) + (pow(g, 2) / (2 * _delta)));
	Operator<T> Ht = pauliZ<T>() * (-delta_q / 2);

	std::complex<double> omega = -2 * (g / _delta)* drive_amplitude;
	//TODO: For Y rotation gate, use Operator<T> Hd = ((a_dagger - a) * i * drive_amplitude) + (pauliY<T>() * omega * 1/2);
	Operator<T> Hd = ((a_dagger + a) * drive_amplitude) + (pauliX<T>() * omega * 1/2);

	Operator<T> H = Hc + Hz + Ht + Hd;

	return H;
}

template<typename T>
std::pair<Operator<T>, double> Transmon<T>::get_Swap_gate(Transmon<T> target, double drive_amplitude, double Wd) {
	double delta_r = abs(Wr - Wd);
	//drive_amplitude = 3.096 * delta_r/(2*PI);

	double delta_a_1 = abs(w - Wd);
	double delta_a_2 = abs(target.w - Wd);

	//Operator<T> a = annihilation(this->photon_number);
	//Operator<T> a_dagger = a.dagger();

	double omega_r_1 = 2 * 2 * PI*drive_amplitude * (this->g / delta_r);
	double omega_r_2 = 2 * 2 * PI*drive_amplitude * (target.g / delta_r);

	double Wa_prime_1 = w + (pow(omega_r_1, 2) / (2 * delta_a_1));
	double Wa_prime_2 = target.w + (pow(omega_r_2, 2) / (2 * delta_a_2));

	double delta_prime_1 = Wa_prime_1 - Wr;
	double delta_prime_2 = Wa_prime_2 - Wr;

	unsigned number_operator = 1;  //a_dagger.transpose().inner(a)
	double Wa_second_1 = Wa_prime_1 + ((2 * (pow(this->g, 2)) / delta_prime_1) * (number_operator + (1 / 2)));
	double Wa_second_2 = Wa_prime_2 + ((2 * (pow(target.g, 2)) / delta_prime_2) * (number_operator + (1 / 2)));
	//double Wa_second_2 = Wa_second_1;
	Operator<T> Hc = pauliI<T>() * Hcavity(delta_r);
	Operator<T> Hz1 = pauliZ<T>() * pauliI<T>() * (Wa_second_1 / 2);
	Operator<T> Hz2 = pauliI<T>() * pauliZ<T>() * (Wa_second_2 / 2);

	Operator<T> sigXA = pauliX<T>() * pauliI<T>();
	Operator<T> sigXB = pauliI<T>() * pauliX<T>();

	Operator<T> sigYA = pauliY<T>() * pauliI<T>();
	Operator<T> sigYB = pauliI<T>() * pauliY<T>();

	std::complex<double> i(0, 1);
	Operator<T> sigYAbackup = sigYA;
	Operator<T> sigYBbackup = sigYB;

	Operator<T> tempA = sigYA * i;
	Operator<T> tempB = sigYB * i;
	Operator<T> sigPlusA = (sigXA + tempA); //pauliPlus<T>() * pauliI<T>();
	Operator<T> sigPlusB = (sigXB + tempB); // pauliI<T>() * pauliPlus<T>();
	sigPlusA * 0.5;
	sigPlusB * 0.5;

	Operator<T> sigMinusA = sigXA - tempA; // pauliMinus<T>() * pauliI<T>();
	Operator<T> sigMinusB = sigXB - tempB; // pauliI<T>() * pauliMinus<T>(); 
	sigMinusA * 0.5;
	sigMinusB * 0.5;

	double coupling = (this->g * target.g * (delta_prime_1 + delta_prime_2)) / (2 * delta_prime_1 * delta_prime_2);

	std::complex<double> J(coupling, 0);
	Operator<T> Hint = (sigPlusA.inner(sigMinusB) + sigMinusA.inner(sigPlusB)) * J;

	Operator<T> H = Hc + Hz1 + Hz2 + Hint;

	double gate_time = (PI / 4) / coupling;

	pair<Operator<T>, double> result(H, gate_time);

	return result;
}
//qubit approximation: two atom levels qubit, j=1, n=1
// t = PI * _delta/(2 * g * drive_amplitude)
// drive_amplitude = 0.03*Wr
template<typename T>
std::pair<Operator<T>, double> Transmon<T>::get_X_Y_gate(double drive_x_amplitude, double drive_y_amplitude, double Wd, double x_theta, double y_theta) {

	double _delta = abs(Wr - this->w);// atom-cavity detuning
	//double n_bar = pow(drive_amplitude, 2) / pow(_delta, 2); //average photon number can be estimated as 0.1	
	//if drive_amplitude not given use estimation n_bar = 0.1
	//double n_bar = 0.1;

	//double delta_r = abs(this->Wr - Wd);
	//double delta_q = abs(this->w - Wd);

	//Operator<T> a = annihilation(this->photon_number);
	//Operator<T> a_dagger = a.dagger();

	//Operator<T> Hc = Hcavity(delta_r);
	//Operator<T> Hz = pauliZ<T>() * (pow(g, 2) / (2 * _delta));
	//Operator<T> Ht = pauliZ<T>() * (-delta_q / 2);

	std::complex<double> omega_x = -2 * (g / _delta)* drive_x_amplitude* (x_theta / PI);
	Operator<T> Hdx = pauliX<T>() * omega_x / 2;

	std::complex<double> omega_y = -2 * (g / _delta)* drive_y_amplitude* (y_theta/PI);
	Operator<T> Hdy = pauliY<T>() * omega_y / 2;

	//Operator<T> H = Hc + Hz + Ht + Hd;
	Operator<T> H = Hdy + Hdx;
	double gate_time = 20;
	//double gate_time = round(theta * _delta / (2 * g * drive_amplitude * 2 * PI));

	pair<Operator<T>, double> result(H, gate_time);

	return result;
}

template<typename T>
Operator<T> Transmon<T>::get_X_or_Y_gate(double drive_amplitude, double Wd, double x_theta, double y_theta) 
{
	double _delta = abs(Wr - this->w);// atom-cavity detuning
	//double n_bar = pow(drive_amplitude, 2) / pow(_delta, 2); //average photon number can be estimated as 0.1	
	//if drive_amplitude not given use estimation n_bar = 0.1
	//double n_bar = 0.1;

	double delta_r = abs(this->Wr - Wd);
	//double delta_q = abs(this->w - Wd);

	//Operator<T> a = annihilation(this->photon_number);
	//Operator<T> a_dagger = a.dagger();

	//Operator<T> Hc = Hcavity(delta_r);
	//Operator<T> Hz = pauliZ<T>() * (pow(g, 2) / (2 * _delta));
	//Operator<T> Ht = pauliZ<T>() * (-delta_q / 2);

	std::complex<double> omega_x = -2 * (g / delta_r)* drive_amplitude ;
	Operator<T> Hdx = pauliX<T>() * x_theta * omega_x/2;

	std::complex<double> omega_y = -2 * (g / delta_r)* drive_amplitude;
	Operator<T> Hdy = pauliY<T>() * y_theta * omega_y/2;

	//Operator<T> H = Hc + Hz + Ht + Hd;
	Operator<T> H = Hdy + Hdx;

	//double gate_time = round(theta * _delta / (2 * g * drive_amplitude * 2 * PI));

	return H;
}

template<typename T>
Operator<T> Transmon<T>::get_Cphase_gate_freq_pulse(Transmon<T> target, double w_l, double w_r) 
{
	this->w = w_l;
	target.w = w_r;

	Operator<T> a_init = annihilation(this->photon_number);
	Operator<T> a_dagger_init = creation(this->photon_number);
	Operator<T> a = pauliI(this->atom_level) * pauliI(this->atom_level) * a_init;//18 by 18
	Operator<T> a_dagger = pauliI(this->atom_level) * pauliI(this->atom_level) * a_dagger_init;//18 by 18

	Operator<T> Hcavity = a_dagger.inner(a) * this->Wr;
	//Hcavity.print();

	Operator<T> H_bare_l_init = get_trans_freq_operator(this->atom_level, w_l);
	//H_bare_l_init.print();

	Operator<T> H_bare_l = H_bare_l_init * pauliI(this->atom_level) * pauliI(this->photon_number);//18 by 18

	Operator<T> H_bare_r_init = get_trans_freq_operator(this->atom_level, w_r);
	Operator<T> H_bare_r = pauliI(this->atom_level) * H_bare_r_init * pauliI(this->photon_number);//18 by 18

	//coupling terms
	Operator<T> c_init = annihilation(this->atom_level);
	//c_init.print();

	Operator<T> cl = c_init * pauliI(this->atom_level) * pauliI(this->photon_number);//18 by 18
	Operator<T> cl_dagger = cl.dagger();//18 by 18
	//cl.print();

	Operator<T> cr = pauliI(this->atom_level) * c_init * pauliI(this->photon_number);//18 by 18
	Operator<T> cr_dagger = cr.dagger();//18 by 18

	Operator<T> coupling_l_init = a_dagger.inner(cl) + a.inner(cl_dagger);
	//coupling_l_init.print();

	Operator<T> coupling_l = coupling_l_init * this->g;
	//coupling_l.print();

	Operator<T> coupling_r_init = a_dagger.inner(cr) + a.inner(cr_dagger);
	Operator<T> coupling_r = coupling_r_init * target.g;
	//coupling_r.print();
		
	////Effective coupling terms
	//double J00 = get_J(0, 0, target.Wr, w_l, w_r, this->delta, target.delta, this->g, target.g);
	//double J01 = get_J(0, 1, target.Wr, w_l, w_r, this->delta, target.delta, this->g, target.g);
	//double J10 = get_J(1, 0, target.Wr, w_l, w_r, this->delta, target.delta, this->g, target.g);
	//double J11 = get_J(1, 1, target.Wr, w_l, w_r, this->delta, target.delta, this->g, target.g);

	//std::deque<Projector<T>> coupling_projectors;

	//Projector<unsigned> proj0({ 0, 1, 0 }, { 1, 0, 0 }, J00);
	//Projector<unsigned> proj1({ 1, 0, 0 }, { 0, 1, 0 }, J00);

	//Projector<unsigned> proj2({ 0, 2, 0 }, { 1, 1, 0 }, sqrt(2) * J01);
	//Projector<unsigned> proj3({ 1, 1, 0 }, { 0, 2, 0 }, sqrt(2) * J01);
	//Projector<unsigned> proj4({ 2, 0, 0 }, { 1, 1, 0 }, sqrt(2) * J10);
	//Projector<unsigned> proj5({ 1, 1, 0 }, { 2, 0, 0 }, sqrt(2) * J10);

	//Projector<unsigned> proj6({ 1, 2, 0 }, { 2, 1, 0 }, 2 * J11);
	//Projector<unsigned> proj7({ 2, 1, 0 }, { 1, 2, 0 }, 2 * J11);

	////Projector<unsigned> proj0(Bra<unsigned>(3, { 0, 1, 0 }, weight), Ket<unsigned>(3, { 1, 0, 0 }, weight), this->g);
	////Projector<unsigned> proj1(Bra<unsigned>(3, { 1, 0, 0 }, weight), Ket<unsigned>(3, { 0, 1, 0 }, weight), this->g);

	////Projector<unsigned> proj2(Bra<unsigned>(3, { 0, 2, 0 }, weight), Ket<unsigned>(3, { 1, 1, 0 }, weight), sqrt(2) * this->g);
	////Projector<unsigned> proj3(Bra<unsigned>(3, { 1, 1, 0 }, weight), Ket<unsigned>(3, {0, 2, 0}, weight), sqrt(2) * this->g);
	////Projector<unsigned> proj4(Bra<unsigned>(3, { 2, 0, 0 }, weight), Ket<unsigned>(3, { 1, 1, 0 }, weight), sqrt(2) * this->g);
	////Projector<unsigned> proj5(Bra<unsigned>(3, { 1, 1, 0 }, weight), Ket<unsigned>(3, { 2, 0, 0 }, weight), sqrt(2) * this->g);

	////Projector<unsigned> proj6(Bra<unsigned>(3, { 1, 2, 0 }, weight), Ket<unsigned>(3, { 2, 1, 0 }, weight), 2 * this->g);
	////Projector<unsigned> proj7(Bra<unsigned>(3, { 2, 1, 0 }, weight), Ket<unsigned>(3, { 1, 2, 0 }, weight), 2 * this->g);

	//coupling_projectors.push_back(proj0);
	//coupling_projectors.push_back(proj1);
	//coupling_projectors.push_back(proj2);
	//coupling_projectors.push_back(proj3);
	//coupling_projectors.push_back(proj4);
	//coupling_projectors.push_back(proj5);
	//coupling_projectors.push_back(proj6);
	//coupling_projectors.push_back(proj7);

	//Operator<T> H_coupling(Hcavity.size, coupling_projectors);

	Operator<T> H_bare = Hcavity + H_bare_l + H_bare_r;
	//H_bare.print();

	Operator<T> H_coupling = coupling_l + coupling_r;
	//H_coupling.print();

	Operator<T> H = H_bare + H_coupling;
	//H.print();

	return H;
}

#endif