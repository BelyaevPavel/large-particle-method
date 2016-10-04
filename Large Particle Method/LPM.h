#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

using namespace std;
class LPM
{
	//Gas and particles constant parameters
	double R; //Gas constant
	double cp; //Gas heat capacity at constant pressure
	double cv; //Gas heat capacity at constant volume
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
	double gasLeftCellEnergyFlow;
	double gasLeftCellVelocityFlow;
	double gasThisCellEnergyFlow;
	double gasThisCellVelocityFlow;
	double particleLeftCellEnergyFlow;
	double particleLeftCellVelocityFlow;
	double particleThisCellEnergyFlow;
	double particleThisCellVelocityFlow;
	double Pr; //Prandtle value
	//Grid parameters
	int N; //Grid cells count
	double Co; //Courant 
	double dx; //Grid step
	double dt; //Time step
	//Arrays declaration
	double *envVelocity;
	double *U; //Solid phase (particles) velocity on current step
	double *UIntermediate; //Solid phase (particles) velocity between steps
	double *U1; //Solid phase (particles) velocity on next step
	double *V; //Gas phase velocity on current step
	double *VIntermediate; //Gas phase velocity between steps
	double *V1; //Gas phase velocity on next step
	double *envE;
	double *gasE; //Gas phase energy on current step
	double *gasEIntermediate; //Gas phase energy between steps
	double *gasE1; //Gas phase energy on next step
	double *particleE; //Gas phase energy on current step
	double *particleEIntermediate; //Gas phase energy between steps
	double *particleE1; //Gas phase energy on next step
	double *envPressure;
	double *gasP; //Gas phase (particles) pressure on current step
	double *gasPIntermediate; //Gas phase (particles) pressure between steps
	double *gasP1; //Gas phase (particles) pressure on next step
	double *particleP; //Gas phase (particles) pressure on current step
	double *particlePIntermediate; //Gas phase (particles) pressure between steps
	double *particleP1; //Gas phase (particles) pressure on next step
	double *envRo;
	double *gasRo; //Gas phase density on current step
	double *gasRo1; //Gas phase density on next step
	double *particleRo; //Gas phase density on current step
	double *particleRo1; //Gas phase density on next step
	double *mixtureS;
	double *gasS;
	double *particleS;
	//Auxiliary arrays declaration
	double *envR;
	double *gasR;
	double *particleR;
	double *envQ;
	double *gasQ;
	double *particleQ;
	double *envF;
	double *gasF;
	double *particleF;
	double *envFi;
	double *gasFi;
	double *particleFi;
	double *alpha; //Gas proportion by volume
	double *eta; //Gas proportion by mass
	double *Re; //Reynolds value
	double *Nu; //Nuselt value
	double *Cd; //Particle's dinamyc viscousity
	double *f; //Interfactional interaction force
	double *q; //Interfactional heat interaction
	double *p_count; //Particles count
	double *p_count1;
	double *viscosity1; //Artificial viscosity
	double *viscosity2; //Artificial viscosity
	double *gasRightBorderMassFlow;
	double *particleRightBorderMassFlow;
	double *RightBorderParticleFlow;
public:
	double DefineFlowDirection(int index, double* IntermVelocity, int Border, double* Value);
	double DefineRightBorderMassTranslation(double ThisCellVelocity, double RightCellVelocity, double ThisCellDensity, double RightCellDensity, double dt);
	void MainProc();
	LPM();
	~LPM();
};

