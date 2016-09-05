#include "LPM.h"


LPM::LPM()
{
	//Gas and particles constant parameters initialisation
	R = 287;
	cp = 1005;
	c2 = 710;
	ro01 = 1.21;
	ro02 = 7800;
	gamma = 1.4;
	lambda1 = 0.026;
	p0 = 0.1e6;
	a10 = 341;
	M = 4.2;
	mu1 = 1.85e-5;
	T0 = 293;
	//Shielding layer constant parameters initialisation
	alpha10 = 0.999;
	d = 0.001;
	lf = 0.45;
	ls = 3;
	lw = 3.5;
	//Gas calculated parameters initialisation
	u1f = 2 / (gamma + 1)*(M - 1 / M)*a10;
	p1 = (1 + 2 * gamma / (gamma + 1)*M*M - 1)*p0;
	ro01f = (gamma + 1)*M*M / (2 + (gamma - 1)*M*M)*ro01;
	//Auxiliary variables initialisation
	a = 0.9;
	b = 0.5;
	xi;
	loopCounter = 0;
	t = 0;
	tstop = 0.01;
	//Grid parameters initialisation
	N = 400;
	Co = 0.1;
	dx = lw / N;
	dt = Co*dx / (u1f*u1f + a10*a10);
	//Arrays declaration initialisation
	U = new double[N];
	UIntermediate = new double[N];
	U1 = new double[N];
	V = new double[N];
	VIntermediate = new double[N];
	V1 = new double[N];
	E = new double[N];
	EIntermediate = new double[N];
	E1 = new double[N];
	P = new double[N];
	PIntermediate = new double[N];
	P1 = new double[N];
	Ro = new double[N];
	Ro1 = new double[N];
	//Auxiliary arrays declaration initialisation
	viscosity1 = new double[N];
	viscosity2 = new double[N];

	for (int i = 0; i < N; i++)
	{
		if (i*dx < lf)
		{
			U[i] = u1f / (i*dx);
			xi = pow(1 - (gamma - 1)*(u1f - U[i]) / (2 * a10), 2 / (gamma - 1));
			P[i] = p1*pow(xi, gamma);
			Ro[i] = ro01f*xi;
			E[i] = P[i] / ((gamma - 1)*Ro[i]) + U[i] * U[i] / 2.0;
		}
		else
		{
			U[i] = 0;
			P[i] = 100000;
			Ro[i] = ro01;
			E[i] = P[i] / ((gamma - 1)*Ro[i]);
		}
	}
}


LPM::~LPM()
{
}


double LPM::DefineFlowDirection(int index, double* IntermVelocity, int Border, double* Value)
{
	if (IntermVelocity[index - 1 + Border] + IntermVelocity[index + Border] == 0)
		return 0;
	else
	{
		if (IntermVelocity[index - 1 + Border] + IntermVelocity[index + Border] > 0)
			return Value[index];
		else
			return Value[index + 1];
	}
}

double LPM::DefineMassTranslatio(double ThisCellVelocity, double RightCellVelocity, double ThisCellDensity, double RightCellDensity, double dt)
{
	if (ThisCellVelocity + RightCellVelocity > 0)
	{
		return ThisCellDensity*(ThisCellDensity + RightCellVelocity)*dt / 2.0;
	}
	else
	{
		if (ThisCellVelocity + RightCellVelocity < 0)
		{
			return RightCellDensity*(ThisCellVelocity + RightCellVelocity)*dt / 2.0;
		}
		else
		{
			return 0;
		}
	}
	return 0;
}

void LPM::MainProc()
{
	// Artificial viscosity calculation.
	for (int i = 2; i < N - 1;i++)
	{
		if (V[i + 1]<V[i])
		{
			viscosity1[i] = ((Ro[i + 1] + Ro[i]) / 2)*(a*sqrt(1.4*(P[i + 1] + P[i]) / (Ro[i + 1] + Ro[i])) + b*(V[i] - V[i + 1]))*(V[i] - V[i + 1]);
			viscosity1[i] = 0;
		}
		if (V[i - 1]>V[i])
		{
			viscosity2[i] = -((Ro[i - 1] + Ro[i]) / 2)*(a*sqrt(1.4*(P[i - 1] + P[i]) / (Ro[i - 1] + Ro[i])) - b*(V[i] - V[i - 1]))*(V[i] - V[i - 1]);
		}
		else
			viscosity2[i] = 0;
	}
	// Intermediate pressure calculation.
	for (int i = 0; i < N-1; i++)
		PIntermediate[i] = (P[i + 1] + P[i]) / 2.0;
	PIntermediate[N - 1] = PIntermediate[N - 2];
	// Intermediate velocity calculation.
	for (int i = 2; i < N - 1; i++)
		VIntermediate[i] = V[i] - ((P[i + 1] - P[i - 1]) / 2.0 + viscosity1[i] - viscosity2[i])*(dt / (dx*Ro[i]));
	VIntermediate[0] = VIntermediate[1];
	VIntermediate[N - 1] = VIntermediate[N - 2];

}