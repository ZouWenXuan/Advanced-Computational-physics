#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>
#include <random>
#include<string>
using namespace std;

// set random seed
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> rdistr(0, 1);

// Model parameters
double J1;			            // ferromagnetic coupling 1
double J2;						// ferromagnetic coupling 2
double T;                       // temperature
int N_spin;                     // number of lattice sites
int* S;                         // the spins

// Monte Carlo parameters
long long N_test;               // number of repeated test
long long N_warmup;             // number of warming up steps
long long N_interval;           // number of MC steps between samples
long long N_samples;            // number of samples of one test

// Observables for one test
double magSum;                  // accumulate M
double m1Sum;					// accumulate M1
double m2Sum;					// accumulate M2
double chiSum;                  // accumulate M^2
double engSum;                  // accumulate E
double engSqSum;                // accumulate E^2
double nAvg;                    // number of values accumulated
double currentE;                // energy of current step

// Observables stored for different test
double* Engy;                   // store energy in each test;
double* Cv;                     // store specific heat in test;
double* Mz;                     // store |Mz| in each test;
double* Mz1;                    // store Mz1 in each test;
double* Mz2;                    // store Mz2 in each test;
double* Xu;                     // store uniform susceptibility in each test;

// time cost
clock_t time1;
clock_t time2;


void initializeSpin() {
	// initialization for the whole test
	// All tests are done on one chain

	// initialize Spin
	S = new int[N_spin];
	for (int i = 0; i < N_spin; i++)
		S[i] = rdistr(gen) < 0.5 ? +1 : -1;

	// initialize current energy 
	double eng = 0.0;
	for (int i = 0; i < N_spin; i++) {
		int iPrev1 = (i - 1 + N_spin) % N_spin;
		int iPrev2 = (i - 2 + N_spin) % N_spin;
		eng += J1 * S[i] * S[iPrev1] + J2 * S[i] * S[iPrev2];
	}
	currentE = eng;

	// initialize observables for all tests
	Engy = new double[N_test];
	Cv = new double[N_test];
	Mz = new double[N_test];
	Mz1 = new double[N_test];
	Mz2 = new double[N_test];
	Xu = new double[N_test];
	for (int i = 0; i < N_test; i++) {
		Engy[i] = 0;
		Cv[i] = 0;
		Mz[i] = 0;
		Mz1[i] = 0;
		Mz2[i] = 0;
		Xu[i] = 0;
	}
}

void initializeObservables() {

	// initialize the sum of obervables for one test
	magSum = 0;
	m1Sum = 0;
	m2Sum = 0;
	chiSum = 0;
	engSum = 0;
	engSqSum = 0;
	nAvg = 0;

}


void singleSpinUpdate() {

	// choose a random spin
	int i = int(rdistr(gen) * N_spin);

	// compute delta energy
	bool accept = false;
	double localfield = J1 * (S[(i - 1 + N_spin) % N_spin] + S[(i + 1 + N_spin) % N_spin]) + J2 * (S[(i - 2 + N_spin) % N_spin] + S[(i + 2 + N_spin) % N_spin]);
	double deltaE = -2 * localfield * S[i];

	// update (or not)
	double prob = exp(-1 / T * deltaE);
	if (deltaE < 0)
	{
		accept = true;
	}
	else {
		if (rdistr(gen) < prob) {
			accept = true;
		}
	}
	if (accept) {
		S[i] = -S[i];
		currentE += deltaE;
	}
}

void measureObservables() {
	// measure observables for the sample in one test

	// observables derived from the magnetic moment
	int M = 0;
	int M1 = 0;
	int M2_odd = 0;
	int M2_even = 0;
	for (int i = 0; i < N_spin; i++) {
		M += S[i];
		M1 += pow(-1, i) * S[i];
		if (i % 2 == 0) {
			M2_even += pow(-1, i / 2) * S[i];
		}
		else {
			M2_odd += pow(-1, (i + 1) / 2) * S[i];
		}

	}

	magSum += double(abs(M));
	chiSum += M * double(M);
	m1Sum += double(abs(M1));
	m2Sum += double(abs(M2_odd));
	m2Sum += double(abs(M2_even));

	// observables derived from the energy moment
	double eng = currentE;
	engSum += eng;
	engSqSum += eng * eng;

	nAvg += 1;

}



void storeObservables(const int i) {
	// store observables of one test

	// For Total spins
	Engy[i] = engSum / double(nAvg);
	Mz[i] = magSum / double(nAvg);
	Mz1[i] = m1Sum / double(nAvg);
	Mz2[i] = m2Sum / double(nAvg);
	Cv[i] = (engSqSum / double(nAvg) - Engy[i] * Engy[i]) / (double(N_spin) * T * T);
	Xu[i] = (chiSum / double(nAvg) - Mz[i] * Mz[i]) / (double(N_spin) * T);

	// Per spins
	Engy[i] /= double(N_spin);
	Mz[i] /= double(N_spin);
	Mz1[i] /= double(N_spin);
	Mz2[i] /= double(N_spin);

}

void computeAverages(string filename) {
	// compute averages for different tests and save

	int i;
	double avgEngy = 0, errEngy = 0;
	double avgMz = 0, errMz = 0;
	double avgMz1 = 0, errMz1 = 0;
	double avgMz2 = 0, errMz2 = 0;
	double avgCv = 0, errCv = 0;
	double avgXu = 0, errXu = 0;
	FILE* fp;

	avgEngy = 0.0;
	for (i = 0; i < N_test; i++)
		avgEngy += Engy[i];
	avgEngy /= double(N_test);

	errEngy = 0;
	for (int i = 0; i < N_test; i++)
		errEngy += (avgEngy - Engy[i]) * (avgEngy - Engy[i]);
	errEngy = pow(errEngy / double(N_test * (N_test)), 0.5);

	avgMz = 0.0;
	for (i = 0; i < N_test; i++)
		avgMz += Mz[i];
	avgMz /= double(N_test);

	errMz = 0;
	for (int i = 0; i < N_test; i++)
		errMz += (avgMz - Mz[i]) * (avgMz - Mz[i]);
	errMz = pow(errMz / double(N_test * (N_test)), 0.5);

	avgMz1 = 0.0;
	for (i = 0; i < N_test; i++)
		avgMz1 += Mz1[i];
	avgMz1 /= double(N_test);

	errMz1 = 0;
	for (int i = 0; i < N_test; i++)
		errMz1 += (avgMz1 - Mz1[i]) * (avgMz1 - Mz1[i]);
	errMz1 = pow(errMz1 / double(N_test * (N_test)), 0.5);

	avgMz2 = 0.0;
	for (i = 0; i < N_test; i++)
		avgMz2 += Mz2[i];
	avgMz2 /= double(N_test);

	errMz2 = 0;
	for (int i = 0; i < N_test; i++)
		errMz2 += (avgMz2 - Mz2[i]) * (avgMz2 - Mz2[i]);
	errMz2 = pow(errMz2 / double(N_test * (N_test)), 0.5);

	avgCv = 0.0;
	for (i = 0; i < N_test; i++)
		avgCv += Cv[i];
	avgCv /= double(N_test);

	errCv = 0;
	for (int i = 0; i < N_test; i++)
		errCv += (avgCv - Cv[i]) * (avgCv - Cv[i]);
	errCv = pow(errCv / double(N_test * (N_test)), 0.5);

	avgXu = 0.0;
	for (i = 0; i < N_test; i++)
		avgXu += Xu[i];
	avgXu /= double(N_test);

	errXu = 0;
	for (int i = 0; i < N_test; i++)
		errXu += (avgXu - Xu[i]) * (avgXu - Xu[i]);
	errXu = pow(errXu / double(N_test * (N_test)), 0.5);


	char name[100];
	strcpy_s(name, filename.c_str());
	fopen_s(&fp, name, "at");
	fprintf(fp, "%24.16e\t%24.16e\t%24.16e\t%24.16e\t%24.16e\t%24.16e\t%24.16e\t%24.16e\t%24.16e\t%24.16e\t%24.16e\t%24.16e\t%24.16e\n", T, avgEngy, errEngy, avgMz, errMz, avgMz1, errMz1, avgMz2, errMz2, avgCv, errCv, avgXu, errXu);
	fclose(fp);

}

void freeMemory() {

	delete[]S;
	delete[]Engy;
	delete[]Mz;
	delete[]Mz1;
	delete[]Mz2;
	delete[]Cv;
	delete[]Xu;

}


int main() {
	
	// model parameter
	N_spin = 20;
	J1 = 1;

	// mc parameter 
	N_test = 5;
	N_warmup = 10000;
	N_interval = 50;
	N_samples = 10000;

	// main loop
	// loop for J2
	for (int cp2 = 17; cp2 < 60; cp2++) {
		J2 = 0.025 * cp2;
		string J2str = to_string(J2);
		J2str = J2str.substr(0, 5);
		string filename = "result_" + J2str + ".txt";

		time1 = clock();
		std::cout << "===========================================================================" << endl;
		std::cout << "Monte Carlo for Anti-Ferromagnetic Ising model, J_1 = 1.000, J_2 = " << J2str << endl;
		std::cout << "===========================================================================" << endl;

		// loop for T
		for (int k = 15; k < 30; k++) {
			
			T = 0.1 + 0.1 * k;
			std::cout << "(" << k + 1 << ")" << " Temperature: " << to_string(T) << ":" << endl;

			// initialize for this T
			initializeSpin();

			// warm up steps
			for (int i = 0; i < N_warmup * N_spin; i++)
				singleSpinUpdate();
			std::cout << "    Finish Warm Up, Begin Sampling..." << endl;

			// loop for repeated test
			for (int i = 0; i < N_test; i++) {
				initializeObservables();

				for (int j = 0; j < N_interval * N_spin * N_samples; j++) {
					singleSpinUpdate();

					if (j % (N_interval * N_spin) == 0) {
						measureObservables();
					}
				}
				storeObservables(i);
				std::cout << "    Repeated samples: " << to_string(i) << " finish..." << endl;
			}

			// compute Averages of differnt samples
			computeAverages(filename);

			freeMemory();

		}
		time2 = clock();
		double endtime = (double)(time2 - time1) / CLOCKS_PER_SEC;
		std::cout << "Total time:" << endtime << endl;
	}
}
