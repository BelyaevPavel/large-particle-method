#pragma once
#include <math.h>
class LPM
{
	//Gas and particles constant parameters
	double R;
	double cp;
	double c2;
	double ro01;
	double ro02;
	double gamma;
	double lambda1;
	double p0;
	double a10;
	double M ;
	double mu1;
	double T0;
	//Shielding layer constant parameters
	double alpha10;
	double d;
	double lf;
	double ls;
	double lw;
	//Gas calculated parameters
	double u1f;
	double p1;
	double ro01f;
	//Auxiliary variables
	double a;
	double b;
	double xi;
	double loopCounter;
	double t;
	double tstop;
	//Grid parameters
	int N;
	double Co;
	double dx;
	double dt;
	//Arrays declaration
	double *U;
	double *UIntermediate;
	double *U1;
	double *V;
	double *VIntermediate;
	double *V1;
	double *E;
	double *EIntermediate;
	double *E1;
	double *P;
	double *PIntermediate;
	double *P1;
	double *Ro;
	double *Ro1;
	//Auxiliary arrays declaration
	double *viscosity1;
	double *viscosity2;
public:
	double DefineFlowDirection(int index, double* IntermVelocity, int Border, double* Value);
	double DefineMassTranslatio(double ThisCellVelocity, double RightCellVelocity, double ThisCellDensity, double RightCellDensity, double dt);
	void MainProc();
	LPM();
	~LPM();
};

