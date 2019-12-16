#include <math.h> 
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "BELLE_efficiency.h"
#include "constants.h"

// // BELLE_D0ToK0Spipi_efficiency_correction.cc
const double x_min = pow(mKs + mPi, 2);
const double x_max = pow(mD0 - mPi, 2);
const double y_min = pow(mPi + mPi, 2);
const double y_max = pow(mD0 - mKs, 2);

inline double lKallen(double x, double y, double z){ return ( x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z );}

double compute_helicity_angle(std::string axisLabel, double qij, double qjk){
	double qi = 0.;
	double qj = 0.;
	double qk = 0.;
	double Q=mD0*mD0;

	// permutate cyclic
	if (axisLabel == "axis_K0SpiRS_pipi") { // K0Spi(RS)+pipi system, i=K0S, j=pi+, k=pi-
		qi = mKs*mKs;
		qj = mPi*mPi;
		qk = mPi*mPi;
	} else if (axisLabel == "axis_pipi_K0SpiWS") { // pipi+K0Spi(WS) system, i=pi+, j=pi-, k=K0S
		qi = mPi*mPi;
		qj = mPi*mPi;
		qk = mKs*mKs;
	} else if (axisLabel == "axis_K0SpiWS_K0SpiRS") { // K0Spi(WS)+K0Spi(RS) system, i=pi-, j=K0S, k=pi+
		qi = mPi*mPi;
		qj = mKs*mKs;
		qk = mPi*mPi;
	} else {
		std::cout << "Wrong axis chosen in compute_helicity_angle. Exiting()" << std::endl;
	}

	return  ( 2*qij*(qj+qk-qjk)+(qij-qi+qj)*(Q-qij-qk) ) / sqrt( lKallen(qij,qi,qj) * lKallen(Q,qij,qk) );
}

/*
 *
 * parameterized efficiency model using helicity angle and mass fitted using Legendre polynomials and needed functions
 *
 */

double compute_legendre_polynomial(double x, int legendre_order)
{
	double legendre_value = 0.;
	if (legendre_order == 0) {
		legendre_value = 1;
	} else if (legendre_order ==  1) {
		legendre_value = x;
	} else if (legendre_order ==  2) {
		legendre_value = (3*x*x-1)/2;
	} else if (legendre_order ==  3) {
		legendre_value = (5*x*x*x - 3*x)/2;
	} else if (legendre_order ==  4) {
		legendre_value = (35*pow(x,4) - 30*x*x + 3)/8;
	} else if (legendre_order ==  5) {
		legendre_value = (63*pow(x,5) - 70*pow(x,3) +15*x)/8;
	} else if (legendre_order ==  6) {
		legendre_value = (231*pow(x,6) - 315*pow(x,4) + 105*x*x -5)/16;
	} else if (legendre_order ==  7) {
		legendre_value = (429*pow(x,7) - 693*pow(x,5) + 315*pow(x,3) -35*x)/16;
	} else if (legendre_order ==  8) {
		legendre_value = (6435*pow(x,8) - 12012*pow(x,6) + 6930*pow(x,4) - 1260 *x*x + 35)/128;
	} else if (legendre_order ==  9) {
		legendre_value = (12155*pow(x,9) - 25740*pow(x,7) + 18018*pow(x,5) - 4620*pow(x,3) + 315*x)/128;
	} else if (legendre_order == 10) {
		legendre_value = (46189*pow(x,10) - 109395*pow(x,8) + 90090 * pow(x,6) - 30030*pow(x,4) + 3465 * x*x - 63)/256;
	} else {
		std::cerr << "Legendre polynomial for n > 10 not available" << std::endl;
		throw;
	}
	return legendre_value;
}

double fitfunction_polynomial_Fermi(double xx, const realVector& par) {
	double polynom_mean0 = par[0];

	double eff0 = par[1];

	double polynom_value0 = (
			1. +
			par[2] * pow(xx - polynom_mean0, 1) +
			par[3] * pow(xx - polynom_mean0, 2) +
			par[4] * pow(xx - polynom_mean0, 3) +
			par[5] * pow(xx - polynom_mean0, 4) +
			par[6] * pow(xx - polynom_mean0, 5) +
			par[7] * pow(xx - polynom_mean0, 6)
			);

	double E_f = par[8];

	double exp_function = exp( E_f / par[9] - xx / par[9] + pow((E_f - xx), 3) * par[10]);

	double return_value = eff0 * polynom_value0 / (1 + 1. / exp_function );

	return return_value;
}

double fitfunction_Chebyshev(double xx, const realVector &p) {
	double fA = x_min;
	double fB = x_max;

	int order = 5;

	realVector fT(order, 0.);
//	double fT[order];

	double x = (2.0 * xx - fA -fB)/(fB-fA);

	if (order == 1) return p[0];
	if (order == 2) return p[0] + x*p[1];

	// build the polynomials
	fT[0] = 1;
	fT[1] = x;
	for (int i = 1; i< order-1; ++i) {
		fT[i+1] =  2 *x * fT[i] - fT[i-1];
	}
	double sum = p[0]*fT[0];
	for (int i = 1; i<= order; ++i) {
		sum += p[i] * fT[i];
	}
	return sum;
}

double function_dalitz_efficiency_mass_01_helicity_01_12(const realVector& x) {
	double mass     = x[0];
	double costheta = x[1];


	// fit results for the parameterized efficiency model
	realVector par_value_coefficient_c0(12,0.);
//	double par_error_coefficient_c0[12];

	par_value_coefficient_c0[0] = 1.94893;
	par_value_coefficient_c0[1] = 0.154412;
	par_value_coefficient_c0[2] = -0.0252537;
	par_value_coefficient_c0[3] = -0.0292038;
	par_value_coefficient_c0[4] = -0.0139714;
	par_value_coefficient_c0[5] = -0.010103;
	par_value_coefficient_c0[6] = 0;
	par_value_coefficient_c0[7] = 0;
	par_value_coefficient_c0[8] = 2.93677;
	par_value_coefficient_c0[9] = 0.340256;
	par_value_coefficient_c0[10] = 3511.03;
	par_value_coefficient_c0[11] = 0;

//	par_error_coefficient_c0[0] = 0;
//	par_error_coefficient_c0[1] = 0.000112038;
//	par_error_coefficient_c0[2] = 0.00223659;
//	par_error_coefficient_c0[3] = 0.00243009;
//	par_error_coefficient_c0[4] = 0.00416152;
//	par_error_coefficient_c0[5] = 0.00273249;
//	par_error_coefficient_c0[6] = 0;
//	par_error_coefficient_c0[7] = 0;
//	par_error_coefficient_c0[8] = 0.00209886;
//	par_error_coefficient_c0[9] = 0.0921559;
//	par_error_coefficient_c0[10] = 326.04;
//	par_error_coefficient_c0[11] = 0;

	realVector par_value_coefficient_c1(6,0.);
//	double par_error_coefficient_c1[6];

	par_value_coefficient_c1[0] = -0.000709756;
	par_value_coefficient_c1[1] = 0.00164371;
	par_value_coefficient_c1[2] = 0.00381571;
	par_value_coefficient_c1[3] = -0.00238563;
	par_value_coefficient_c1[4] = 0.000924276;
	par_value_coefficient_c1[5] = -0.000100673;

//	par_error_coefficient_c1[0] = 0.00017585;
//	par_error_coefficient_c1[1] = 0.000334433;
//	par_error_coefficient_c1[2] = 0.000282188;
//	par_error_coefficient_c1[3] = 0.000276603;
//	par_error_coefficient_c1[4] = 0.000209574;
//	par_error_coefficient_c1[5] = 0.000205315;

	realVector par_value_coefficient_c2(6,0.);
//	double par_error_coefficient_c2[6];

	par_value_coefficient_c2[0] = -0.000654936;
	par_value_coefficient_c2[1] = -0.00200889;
	par_value_coefficient_c2[2] = -0.0018105;
	par_value_coefficient_c2[3] = -0.000791849;
	par_value_coefficient_c2[4] = -0.00273026;
	par_value_coefficient_c2[5] = 0.00022377;

//	par_error_coefficient_c2[0] = 0.000226047;
//	par_error_coefficient_c2[1] = 0.000429108;
//	par_error_coefficient_c2[2] = 0.000362395;
//	par_error_coefficient_c2[3] = 0.000355053;
//	par_error_coefficient_c2[4] = 0.000269661;
//	par_error_coefficient_c2[5] = 0.000264444;

	realVector par_value_coefficient_c3(6,0.);
//	double par_error_coefficient_c3[6];

	par_value_coefficient_c3[0] = 0.00202458;
	par_value_coefficient_c3[1] = -0.000606997;
	par_value_coefficient_c3[2] = -0.00185459;
	par_value_coefficient_c3[3] = 0.000334898;
	par_value_coefficient_c3[4] = -0.00155119;
	par_value_coefficient_c3[5] = 0.00023199;

//	par_error_coefficient_c3[0] = 0.00026861;
//	par_error_coefficient_c3[1] = 0.000510008;
//	par_error_coefficient_c3[2] = 0.000430774;
//	par_error_coefficient_c3[3] = 0.000422083;
//	par_error_coefficient_c3[4] = 0.000320603;
//	par_error_coefficient_c3[5] = 0.00031426;

	realVector par_value_coefficient_c4(6.,0.);
//	double par_error_coefficient_c4[6];

	par_value_coefficient_c4[0] = -0.00228998;
	par_value_coefficient_c4[1] = -0.001751;
	par_value_coefficient_c4[2] = -0.0024826;
	par_value_coefficient_c4[3] = -5.89708e-05;
	par_value_coefficient_c4[4] = -0.00158534;
	par_value_coefficient_c4[5] = -0.000672005;

//	par_error_coefficient_c4[0] = 0.000305997;
//	par_error_coefficient_c4[1] = 0.000580573;
//	par_error_coefficient_c4[2] = 0.000490491;
//	par_error_coefficient_c4[3] = 0.000480409;
//	par_error_coefficient_c4[4] = 0.000365158;
//	par_error_coefficient_c4[5] = 0.000358084;

	realVector par_value_coefficient_c5(6,0.);
//	double par_error_coefficient_c5[6];

	par_value_coefficient_c5[0] = 0.0023321;
	par_value_coefficient_c5[1] = 0.00181519;
	par_value_coefficient_c5[2] = 0.000192932;
	par_value_coefficient_c5[3] = 0.000562951;
	par_value_coefficient_c5[4] = 0.000261524;
	par_value_coefficient_c5[5] = 0.000100477;

//	par_error_coefficient_c5[0] = 0.000342731;
//	par_error_coefficient_c5[1] = 0.000650813;
//	par_error_coefficient_c5[2] = 0.000549616;
//	par_error_coefficient_c5[3] = 0.000538791;
//	par_error_coefficient_c5[4] = 0.000409173;
//	par_error_coefficient_c5[5] = 0.000401268;

	realVector par_value_coefficient_c6(6,0.);
//	double par_error_coefficient_c6[6];

	par_value_coefficient_c6[0] = -0.00202768;
	par_value_coefficient_c6[1] = -0.000127716;
	par_value_coefficient_c6[2] = -0.00272417;
	par_value_coefficient_c6[3] = 0.000554799;
	par_value_coefficient_c6[4] = -0.00117353;
	par_value_coefficient_c6[5] = 0.000492697;

//	par_error_coefficient_c6[0] = 0.000378692;
//	par_error_coefficient_c6[1] = 0.000719894;
//	par_error_coefficient_c6[2] = 0.000607639;
//	par_error_coefficient_c6[3] = 0.00059541;
//	par_error_coefficient_c6[4] = 0.000451484;
//	par_error_coefficient_c6[5] = 0.000442199;

	realVector par_value_coefficient_c7(6,0.);
//	double par_error_coefficient_c7[6];

	par_value_coefficient_c7[0] = 0.00176227;
	par_value_coefficient_c7[1] = 0.00118788;
	par_value_coefficient_c7[2] = 0.000404963;
	par_value_coefficient_c7[3] = 0.00051246;
	par_value_coefficient_c7[4] = 0.000331268;
	par_value_coefficient_c7[5] = -4.13886e-05;

//	par_error_coefficient_c7[0] = 0.000416703;
//	par_error_coefficient_c7[1] = 0.000792841;
//	par_error_coefficient_c7[2] = 0.000668794;
//	par_error_coefficient_c7[3] = 0.000655354;
//	par_error_coefficient_c7[4] = 0.000496344;
//	par_error_coefficient_c7[5] = 0.000485863;


	// compute the legendre coefficienct as a function of the mass
	double c0 = fitfunction_polynomial_Fermi( mass, par_value_coefficient_c0 );
	double c1 = fitfunction_Chebyshev( mass, par_value_coefficient_c1);
	double c2 = fitfunction_Chebyshev( mass, par_value_coefficient_c2);
	double c3 = fitfunction_Chebyshev( mass, par_value_coefficient_c3);
	double c4 = fitfunction_Chebyshev( mass, par_value_coefficient_c4);
	double c5 = fitfunction_Chebyshev( mass, par_value_coefficient_c5);
	double c6 = fitfunction_Chebyshev( mass, par_value_coefficient_c6);
	double c7 = fitfunction_Chebyshev( mass, par_value_coefficient_c7);
	double efficiency = c0 * compute_legendre_polynomial(costheta, 0)
		+ c1 * compute_legendre_polynomial(costheta, 1)
		+ c2 * compute_legendre_polynomial(costheta, 2)
		+ c3 * compute_legendre_polynomial(costheta, 3)
		+ c4 * compute_legendre_polynomial(costheta, 4)
		+ c5 * compute_legendre_polynomial(costheta, 5)
		+ c6 * compute_legendre_polynomial(costheta, 6)
		+ c7 * compute_legendre_polynomial(costheta, 7);
	return efficiency;
}

double Efficiency(const realVector& x) {

// x[0] is the mass of the decaying particle (Here, D^0 mass)


	double mass_squared01 = x[1];
	double mass_squared12 = x[2];


	double helicity_costheta_01_12 = compute_helicity_angle( "axis_K0SpiRS_pipi", mass_squared01, mass_squared12 );


	// check if in the Dalitz plane
	if (mass_squared01 < pow(mKs + mPi, 2)) return 0;
	if (pow(mD0 - mPi, 2) < mass_squared01) return 0;
	if (helicity_costheta_01_12 < -1) return 0;
	if (helicity_costheta_01_12 > 1) return 0;


	realVector array_mass_costheta(2,0.);
	array_mass_costheta[0] = mass_squared01;
	array_mass_costheta[1] = helicity_costheta_01_12;

//	double array_par_dummy[50];
	double efficiency = function_dalitz_efficiency_mass_01_helicity_01_12(array_mass_costheta);
	return efficiency;
}
