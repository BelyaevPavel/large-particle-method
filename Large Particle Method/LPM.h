#pragma once
#include <math.h>
#include <iostream>

using namespace std;
class LPM
{
	//Gas and particles constant parameters
	double R; //Gas constant
	double cp; //Gas heat capacity at constant pressure
	double c2; //Particles material heat capacity
	double ro01; //Gas density
	double ro02; //Particles density
	double gamma; //Adiabatic exponent
	double lambda1; //Thermal conduction
	double p0; //Pressure stagnation value
	double a10; //Sound speed stagnation value
	double a1; //
	double M; //Mach
	double mu1; //Viscousity
	double T0; //Temperature stagnation value
	//Shielding layer constant parameters
	double alpha10; //Gas content percente
	double d; //Particle's diameter
	double lf;
	double ls;
	double lw;
	//Gas calculated parameters
	double u1f;
	double p1;
	double ro01f;
	double e0;
	//Auxiliary variables
	double a;
	double b;
	double xi;
	int loopCounter;
	double t; //Current time
	double tstop; //Calculation stop time
	double LeftCellEnergyFlow;
	double LeftCellVelocityFlow;
	double ThisCellEnergyFlow;
	double ThisCellVelocityFlow;
	//Grid parameters
	int N; //Grid cells count
	double Co; //Courant 
	double dx; //Grid step
	double dt; //Time step
	//Arrays declaration
	double *U; //Solid phase (particles) velocity on current step
	double *UIntermediate; //Solid phase (particles) velocity between steps
	double *U1; //Solid phase (particles) velocity on next step
	double *V; //Gas phase (particles) velocity on current step
	double *VIntermediate; //Gas phase (particles) velocity between steps
	double *V1; //Gas phase (particles) velocity on next step
	double *E; //Gas phase (particles) energy on current step
	double *EIntermediate; //Gas phase (particles) energy between steps
	double *E1; //Gas phase (particles) energy on next step
	double *P;  //Gas phase (particles) pressure on current step
	double *PIntermediate; //Gas phase (particles) pressure between steps
	double *P1; //Gas phase (particles) pressure on next step
	double *Ro;  //Gas phase (particles) density on current step
	double *Ro1;  //Gas phase (particles) density on next step
	//Auxiliary arrays declaration
	double *viscosity1; //Artificial viscosity
	double *viscosity2; //Artificial viscosity
	double *RightBorderMassFlow;
public:
	double DefineFlowDirection(int index, double* IntermVelocity, int Border, double* Value);
	double DefineRightBorderMassTranslation(double ThisCellVelocity, double RightCellVelocity, double ThisCellDensity, double RightCellDensity, double dt);
	void MainProc();
	LPM();
	~LPM();
};

