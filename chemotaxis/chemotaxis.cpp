/*#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

#include "vector.hpp"
#include "odeint.hpp"

using cpl::Vector;

enum Variables {
	E0U = 1,
	R = 2,
	E0U_R = 3,
	E1U = 4,
	B = 5,
	E1U_B = 6,
	E0O = 7,
	E0O_R = 8,
	E1O = 9,
	E1O_B = 10,
	l = 11,
	E1U_R = 12,
	E2U = 13,
	E2U_B = 14,
	E1O_R = 15,
	E2O = 16,
	E2O_B = 17,
	E2U_R = 18,
	E3U = 19,
	E3U_B = 20,
	E2O_R = 21,
	E3O = 22,
	E3O_B = 23,
	E3U_R = 24,
	E4U = 25,
	E4U_B = 26,
	E3O_R = 27,
	E4O = 28,
	E4O_B = 29
};

enum class SimulationMode {
	SteadyState,
	Constant,
	Gradient
};

SimulationMode mode = SimulationMode::SteadyState;
double c = 0;

double accuracy = 1e-6;
double dataTimeInterval = 10;

double getLigandConcentration (double time)
{
	switch (mode) {
	case SimulationMode::SteadyState:
		return 0;
	case SimulationMode::Constant:
		return c;
	case SimulationMode::Gradient:
		return time / 10000000.;

	}
	
}

Vector get_derivs (const Vector& conc)
{
	Vector derivs (conc.dimension ());

	const double L = getLigandConcentration (conc[0]);

	const double E2O__R = conc[E2O] * conc[R];
	const double E1U__R = conc[E1U] * conc[R];
	const double E2U__B = conc[E2U] * conc[B];
	const double E2U_B__l = conc[E2U_B] * L;
	const double E1U__l = conc[E1U] * L;
	const double E4U__B = conc[E4U] * conc[B];
	const double E2U__R = conc[E2U] * conc[R];
	const double E0U__l = conc[E0U] * L;
	const double E3U__R = conc[E3U] * conc[R];
	const double E1U_B__l = conc[E1U_B] * L;
	const double E3U__l = conc[E3U] * L;
	const double E1O__R = conc[E1O] * conc[R];
	const double E0O__R = conc[E0O] * conc[R];
	const double E3O__B = conc[E3O] * conc[B];
	const double E4O__B = conc[E4O] * conc[B];
	const double E0U_R__l = conc[E0U_R] * L;
	const double E4U__l = conc[E4U] * L;
	const double E0U__R = conc[E0U] * conc[R];
	const double E4U_B__l = conc[E4U_B] * L;
	const double E3O__R = conc[E3O] * conc[R];
	const double E2U_R__l = conc[E2U_R] * L;
	const double E1O__B = conc[E1O] * conc[B];
	const double E1U_R__l = conc[E1U_R] * L;
	const double E3U_B__l = conc[E3U_B] * L;
	const double E2U__l = conc[E2U] * L;
	const double E3U__B = conc[E3U] * conc[B];
	const double E1U__B = conc[E1U] * conc[B];
	const double E2O__B = conc[E2O] * conc[B];
	const double E3U_R__l = conc[E3U_R] * L;

	derivs[0] = 1;//time

	derivs[E0U] = conc[E0U_R] * 0.1 + conc[E1U_B] * 0.0001 + conc[E0O] * 1 - E0U__R * 0.08 - E0U__l * 1;
	derivs[R] = conc[E0U_R] * 0.1 + conc[E0U_R] * 0.0001 + conc[E0O_R] * 0.1 + conc[E0O_R] * 0.0001 + conc[E1U_R] * 0.1 + conc[E1U_R] * 0.0001 + conc[E1O_R] * 0.1 + conc[E1O_R] * 0.0001 + conc[E2U_R] * 0.1 + conc[E2U_R] * 0.0001 + conc[E2O_R] * 0.1 + conc[E2O_R] * 0.0001 + conc[E3U_R] * 0.1 + conc[E3U_R] * 0.0001 + conc[E3O_R] * 0.1 + conc[E3O_R] * 0.0001 - E0U__R * 0.08 - E0O__R * 0.08 - E1U__R * 0.08 - E1O__R * 0.08 - E2U__R * 0.08 - E2O__R * 0.08 - E3U__R * 0.08 - E3O__R * 0.08;
	derivs[E0U_R] = E0U__R * 0.08 + conc[E0O_R] * 1 - conc[E0U_R] * 0.1 - conc[E0U_R] * 0.0001 - E0U_R__l * 1;
	derivs[E1U] = conc[E0U_R] * 0.0001 + conc[E1U_B] * 1.0 + conc[E1U_R] * 0.1 + conc[E2U_B] * 0.0001 + conc[E1O] * 1 - E1U__B * 0.08000000000000002 - E1U__R * 0.08 - E1U__l * 1;
	derivs[B] = conc[E1U_B] * 1.0 + conc[E1U_B] * 0.0001 + conc[E1O_B] * 1.0 + conc[E1O_B] * 0.0001 + conc[E2U_B] * 1.0 + conc[E2U_B] * 0.0001 + conc[E2O_B] * 1.0 + conc[E2O_B] * 0.0001 + conc[E3U_B] * 1.0 + conc[E3U_B] * 0.0001 + conc[E3O_B] * 1.0 + conc[E3O_B] * 0.0001 + conc[E4U_B] * 1.0 + conc[E4U_B] * 0.0001 + conc[E4O_B] * 1.0 + conc[E4O_B] * 0.0001 - E1U__B * 0.08000000000000002 - E1O__B * 0.0 - E2U__B * 0.4 - E2O__B * 0.08000000000000002 - E3U__B * 0.6000000000000001 - E3O__B * 0.4 - E4U__B * 0.8 - E4O__B * 0.8;
	derivs[E1U_B] = E1U__B * 0.08000000000000002 + conc[E1O_B] * 1 - conc[E1U_B] * 1.0 - conc[E1U_B] * 0.0001 - E1U_B__l * 1;
	derivs[E0O] = conc[E0O_R] * 0.1 + conc[E1O_B] * 0.0001 + E0U__l * 1 - E0O__R * 0.08 - conc[E0O] * 1;
	derivs[E0O_R] = E0O__R * 0.08 + E0U_R__l * 1 - conc[E0O_R] * 0.1 - conc[E0O_R] * 0.0001 - conc[E0O_R] * 1;
	derivs[E1O] = conc[E0O_R] * 0.0001 + conc[E1O_B] * 1.0 + conc[E1O_R] * 0.1 + conc[E2O_B] * 0.0001 + E1U__l * 1 - E1O__B * 0.0 - E1O__R * 0.08 - conc[E1O] * 1;
	derivs[E1O_B] = E1O__B * 0.0 + E1U_B__l * 1 - conc[E1O_B] * 1.0 - conc[E1O_B] * 0.0001 - conc[E1O_B] * 1;
	derivs[l] = conc[E0O] * 1 + conc[E0O_R] * 1 + conc[E1O_B] * 1 + conc[E1O] * 1 + conc[E1O_R] * 1 + conc[E2O_B] * 1 + conc[E2O] * 1 + conc[E2O_R] * 1 + conc[E3O_B] * 1 + conc[E3O] * 1 + conc[E3O_R] * 1 + conc[E4O_B] * 1 + conc[E4O] * 1 - E0U__l * 1 - E0U_R__l * 1 - E1U_B__l * 1 - E1U__l * 1 - E1U_R__l * 1 - E2U_B__l * 1 - E2U__l * 1 - E2U_R__l * 1 - E3U_B__l * 1 - E3U__l * 1 - E3U_R__l * 1 - E4U_B__l * 1 - E4U__l * 1;
	derivs[E1U_R] = E1U__R * 0.08 + conc[E1O_R] * 1 - conc[E1U_R] * 0.1 - conc[E1U_R] * 0.0001 - E1U_R__l * 1;
	derivs[E2U] = conc[E1U_R] * 0.0001 + conc[E2U_B] * 1.0 + conc[E2U_R] * 0.1 + conc[E3U_B] * 0.0001 + conc[E2O] * 1 - E2U__B * 0.4 - E2U__R * 0.08 - E2U__l * 1;
	derivs[E2U_B] = E2U__B * 0.4 + conc[E2O_B] * 1 - conc[E2U_B] * 1.0 - conc[E2U_B] * 0.0001 - E2U_B__l * 1;
	derivs[E1O_R] = E1O__R * 0.08 + E1U_R__l * 1 - conc[E1O_R] * 0.1 - conc[E1O_R] * 0.0001 - conc[E1O_R] * 1;
	derivs[E2O] = conc[E1O_R] * 0.0001 + conc[E2O_B] * 1.0 + conc[E2O_R] * 0.1 + conc[E3O_B] * 0.0001 + E2U__l * 1 - E2O__B * 0.08000000000000002 - E2O__R * 0.08 - conc[E2O] * 1;
	derivs[E2O_B] = E2O__B * 0.08000000000000002 + E2U_B__l * 1 - conc[E2O_B] * 1.0 - conc[E2O_B] * 0.0001 - conc[E2O_B] * 1;
	derivs[E2U_R] = E2U__R * 0.08 + conc[E2O_R] * 1 - conc[E2U_R] * 0.1 - conc[E2U_R] * 0.0001 - E2U_R__l * 1;
	derivs[E3U] = conc[E2U_R] * 0.0001 + conc[E3U_B] * 1.0 + conc[E3U_R] * 0.1 + conc[E4U_B] * 0.0001 + conc[E3O] * 1 - E3U__B * 0.6000000000000001 - E3U__R * 0.08 - E3U__l * 1;
	derivs[E3U_B] = E3U__B * 0.6000000000000001 + conc[E3O_B] * 1 - conc[E3U_B] * 1.0 - conc[E3U_B] * 0.0001 - E3U_B__l * 1;
	derivs[E2O_R] = E2O__R * 0.08 + E2U_R__l * 1 - conc[E2O_R] * 0.1 - conc[E2O_R] * 0.0001 - conc[E2O_R] * 1;
	derivs[E3O] = conc[E2O_R] * 0.0001 + conc[E3O_B] * 1.0 + conc[E3O_R] * 0.1 + conc[E4O_B] * 0.0001 + E3U__l * 1 - E3O__B * 0.4 - E3O__R * 0.08 - conc[E3O] * 1;
	derivs[E3O_B] = E3O__B * 0.4 + E3U_B__l * 1 - conc[E3O_B] * 1.0 - conc[E3O_B] * 0.0001 - conc[E3O_B] * 1;
	derivs[E3U_R] = E3U__R * 0.08 + conc[E3O_R] * 1 - conc[E3U_R] * 0.1 - conc[E3U_R] * 0.0001 - E3U_R__l * 1;
	derivs[E4U] = conc[E3U_R] * 0.0001 + conc[E4U_B] * 1.0 + conc[E4O] * 1 - E4U__B * 0.8 - E4U__l * 1;
	derivs[E4U_B] = E4U__B * 0.8 + conc[E4O_B] * 1 - conc[E4U_B] * 1.0 - conc[E4U_B] * 0.0001 - E4U_B__l * 1;
	derivs[E3O_R] = E3O__R * 0.08 + E3U_R__l * 1 - conc[E3O_R] * 0.1 - conc[E3O_R] * 0.0001 - conc[E3O_R] * 1;
	derivs[E4O] = conc[E3O_R] * 0.0001 + conc[E4O_B] * 1.0 + E4U__l * 1 - E4O__B * 0.8 - conc[E4O] * 1;
	derivs[E4O_B] = E4O__B * 0.8 + E4U_B__l * 1 - conc[E4O_B] * 1.0 - conc[E4O_B] * 0.0001 - conc[E4O_B] * 1;

	return derivs;
}


Vector initialValue ()
{
	Vector conc (30);

	conc[0] = 0;

	conc[E0U] = 11.861270579304454;
	conc[R] = 0.23722541158608912;
	conc[E0U_R] = 0;
	conc[E1U] = 0;
	conc[B] = 2.372254115860891;
	conc[E1U_B] = 0;
	conc[E0O] = 0;
	conc[E0O_R] = 0;
	conc[E1O] = 0;
	conc[E1O_B] = 0;
	conc[l] = 0;
	conc[E1U_R] = 0;
	conc[E2U] = 0;
	conc[E2U_B] = 0;
	conc[E1O_R] = 0;
	conc[E2O] = 0;
	conc[E2O_B] = 0;
	conc[E2U_R] = 0;
	conc[E3U] = 0;
	conc[E3U_B] = 0;
	conc[E2O_R] = 0;
	conc[E3O] = 0;
	conc[E3O_B] = 0;
	conc[E3U_R] = 0;
	conc[E4U] = 0;
	conc[E4U_B] = 0;
	conc[E3O_R] = 0;
	conc[E4O] = 0;
	conc[E4O_B] = 0;

	return conc;
}

double getActivity (const Vector& conc)
{
	const double a_m[] = { 0, 0.1, 0.5, 0.75, 1 };
	const double a_m_O[] = { 0, 0, 0.1, 0.5, 1 };
	return (a_m[0] * conc[E0U] + a_m[1] * conc[E1U] + a_m[2] * conc[E2U] + a_m[3] * conc[E3U] + a_m[4] * conc[E4U]) +
		(a_m_O[0] * conc[E0O] + a_m_O[1] * conc[E1O] + a_m_O[2] * conc[E2O] + a_m_O[3] * conc[E3O] + a_m_O[4] * conc[E4O]);
}

Vector runSimulation (const Vector& x0, double time, std::vector<Vector>* data = nullptr)
{
	Vector x = x0;

	const double until = x0[0] + time;
	double dt = 0.01;
	double lastPlotTime = -1;

	double dt_max = 0, dt_min = 100;

	auto start = std::chrono::system_clock::now ();

	int steps = 0;

	std::cout << "Starting new simulation phase" << std::endl;
	do {
		//if (data != nullptr && x[0] > lastPlotTime + dataTimeInterval) {
		if (steps % 100 == 0) {
			lastPlotTime = x[0];
			data->push_back (x);
		}
		
		adaptiveRK4Step (x, dt, accuracy, get_derivs);

		if (steps % 100000 == 0)
			std::cout << "time: " << x[0] - x0[0] << " / " << time << std::endl;

		steps++;
		if (dt < dt_min) dt_min = dt;
		if (dt > dt_max) dt_max = dt;
	} while (x[0] < until);

	auto end = std::chrono::system_clock::now ();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count ();

	std::cout << "Runtime\t" << elapsed << " ms" << std::endl;
	std::cout << " number of adaptive steps = " << steps << std::endl;
	std::cout << " step size: min = " << dt_min << "  max = " << dt_max << std::endl;

	return x;
}

void writeDataToFile (const std::string& filename, const std::vector<Vector>& data)
{
	std::ofstream dataFile (filename);

	const double startTime = data[0][0];

	for (const Vector& x : data) {
		dataFile << x[0] - startTime << "\t" << getActivity (x) << "\t";
		
		for (int i = 1; i < x.dimension (); i++)
			dataFile << x[i] << '\t';

		dataFile << '\n';
	}

	dataFile.close ();
}

void makeSingleStep (const Vector& x0, double L)
{
	std::vector<Vector> data;

	mode = SimulationMode::Constant;
	c = L;
	runSimulation (x0, 1e6, &data);

	writeDataToFile ("chemo_step.dat", data);
}

void makeSteps (const Vector& x0, double c0, const std::vector<double>& ligandConcentrations, double deltaT = 20*60*1000)
{
	std::vector<Vector> data;

	Vector x = x0;

	mode = SimulationMode::Constant;
	c = c0;

	for (double l : ligandConcentrations) {
		c += l;
		x = runSimulation (x, deltaT, &data);
	}

	writeDataToFile ("steps.dat", data);
}

void makeGradient (const Vector& x0)
{
	std::vector<Vector> data;
	
	Vector x = x0;
	x[0] = 0;

	mode = SimulationMode::Gradient;

	runSimulation (x, 1e7, &data);

	writeDataToFile ("gradient.dat", data);
}

int main ()
{
	std::vector<Vector> steady_data;

	const Vector x_steady = runSimulation (initialValue (), 5e5, &steady_data);
	
	writeDataToFile ("chemo_steady.dat", steady_data);

	makeSingleStep (x_steady, 0.1);

	makeSteps (x_steady, 0, {1, -1, 3, -3, 5, -5, 7, -7});

	makeGradient (x_steady);

	system ("PAUSE");
}*/