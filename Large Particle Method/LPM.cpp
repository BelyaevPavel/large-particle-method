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
	p1 = (1 + 2 * gamma*(M*M - 1) / (gamma + 1))*p0;//(1 + 2 * gamma / (gamma + 1)*M*M - 1)*p0;
	ro01f = (gamma + 1)*M*M / (2 + (gamma - 1)*M*M)*ro01;
	a10 = sqrt(gamma*p0 / ro01);
	a1 = sqrt(gamma*p1 / ro01f);
	u1f = 2 / (gamma + 1)*(M - 1 / M)*a10;
	e0 = p0 / (ro01*(gamma - 1));
	//Auxiliary variables initialisation
	a = 0.9;
	b = 0.5;
	xi;
	loopCounter = 0;
	t = 0;
	tstop = 0.01;
	//Grid parameters initialisation
	N = 400;
	Co = 0.3;
	dx = lw / N;
	dt = Co*dx / (sqrt(u1f*u1f + a10*a10));
	//Arrays initialisation
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
	//Auxiliary arrays initialisation
	viscosity1 = new double[N];
	viscosity2 = new double[N];
	RightBorderMassFlow = new double[N];

	for (int i = 0; i < N; i++)
	{
		if (i*dx < lf)
		{
			V[i] = u1f * (i*dx) / lf;
			xi = pow(1 - (gamma - 1)*(u1f - V[i]) / (2 * a10), 2 / (gamma - 1));
			P[i] = p1*pow(xi, gamma);
			Ro[i] = ro01f*xi;
			//V[i] = u1f;
			//P[i] = p1;
			//Ro[i] = ro01f;
			E[i] = P[i] / ((gamma - 1)*Ro[i]) + V[i] * V[i] / 2.0;
		}
		else
		{
			V[i] = 0;
			P[i] = 100000;
			Ro[i] = ro01;
			E[i] = e0;
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

double LPM::DefineRightBorderMassTranslation(double ThisCellVelocity, double RightCellVelocity, double ThisCellDensity, double RightCellDensity, double dt)
{
	if (ThisCellVelocity + RightCellVelocity > 0)
	{
		return ThisCellDensity*(ThisCellVelocity + RightCellVelocity)*dt / 2.0;
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
	double* temp;
	int n = this->N;
	while (t <= tstop)
	{
		loopCounter++;
		// Artificial viscosity calculation.
		for (int i = 1; i < N - 1; i++)
		{
			if (V[i + 1]<V[i])
				viscosity1[i] = ((Ro[i + 1] + Ro[i]) / 2)*(a*sqrt(1.4*(P[i + 1] + P[i]) / (Ro[i + 1] + Ro[i])) + b*(V[i] - V[i + 1]))*(V[i] - V[i + 1]);
			else
				viscosity1[i] = 0;

			if (V[i - 1]>V[i])
				viscosity2[i] = -((Ro[i - 1] + Ro[i]) / 2)*(a*sqrt(1.4*(P[i - 1] + P[i]) / (Ro[i - 1] + Ro[i])) - b*(V[i] - V[i - 1]))*(V[i] - V[i - 1]);
			else
				viscosity2[i] = 0;
		}
		viscosity1[0] = viscosity1[1];
		viscosity1[N - 1] = viscosity1[N - 2];
		viscosity2[0] = viscosity2[1];
		viscosity2[N - 1] = viscosity2[N - 2];
		// Intermediate pressure calculation.
		for (int i = 0; i < N - 1; i++)
		{
			PIntermediate[i] = (P[i + 1] + P[i]) / 2.0;
			if (PIntermediate[i] < 0)
				cout << loopCounter << "Fail." << endl;
		}
		PIntermediate[N - 1] = PIntermediate[N - 2];
		// Intermediate velocity calculation.
		for (int i = 1; i < N - 1; i++)
		{
			VIntermediate[i] = V[i] - ((P[i + 1] - P[i - 1]) / 2.0 + viscosity1[i] - viscosity2[i])*(dt / (dx*Ro[i]));

			//if (VIntermediate[i] < 0)
			//	cout << loopCounter << "Fail." << endl;
		}
		VIntermediate[0] = VIntermediate[1];
		VIntermediate[N - 1] = -VIntermediate[N - 2];
		// Intermediate energy calculation
		for (int i = 1; i < N - 1; i++)
		{
			EIntermediate[i] = E[i] - ((PIntermediate[i] + viscosity1[i])*(V[i] + V[i + 1]) / 2 - (PIntermediate[i - 1] + viscosity2[i])*(V[i] + V[i - 1]) / 2)*(dt / (dx*Ro[i]));

			if (EIntermediate[i] < 0)
				cout << loopCounter << "Fail." << endl;
		}
		EIntermediate[0] = EIntermediate[1];
		EIntermediate[N - 1] = EIntermediate[N - 2];
		// Calculation of mass flow through right cell border
		for (int i = 0; i < N - 1; i++)
			RightBorderMassFlow[i] = DefineRightBorderMassTranslation(V[i], V[i + 1], Ro[i], Ro[i + 1], dt);
		RightBorderMassFlow[N] = 0;
		// Next time step density calculation
		for (int i = 1; i < N; i++)
			Ro1[i] = Ro[i] + (RightBorderMassFlow[i - 1] - RightBorderMassFlow[i]) / dx;
		Ro1[0] = Ro1[1];
		// Final calculations cycle
		for (int i = 1; i < N - 1; i++)
		{
			LeftCellVelocityFlow = DefineFlowDirection(i - 1, VIntermediate, 1, VIntermediate);
			LeftCellEnergyFlow = DefineFlowDirection(i - 1, VIntermediate, 1, EIntermediate);
			ThisCellEnergyFlow = DefineFlowDirection(i, VIntermediate, 1, EIntermediate);
			ThisCellVelocityFlow = DefineFlowDirection(i, VIntermediate, 1, VIntermediate);
			E1[i] = (Ro[i] * dx) / (Ro1[i] * dx)*EIntermediate[i] + (1.0 / (Ro1[i] * dx))*(LeftCellEnergyFlow * RightBorderMassFlow[i - 1] - ThisCellEnergyFlow * RightBorderMassFlow[i]);
			V1[i] = (Ro[i] * dx) / (Ro1[i] * dx)*VIntermediate[i] + (1.0 / (Ro1[i] * dx))*(LeftCellVelocityFlow * RightBorderMassFlow[i - 1] - ThisCellVelocityFlow * RightBorderMassFlow[i]);
			P1[i] = Ro1[i] * (E1[i] - V1[i] * V1[i] / 2.0)*0.4;
			//if (V1[i] < 0)
			//	cout << loopCounter << "Fail." << endl;
		}
		E1[0] = E1[1];
		E1[N - 1] = E1[N - 2];
		P1[0] = P1[1];
		P1[N - 1] = P1[N - 2];
		V1[0] = V1[1];
		V1[N - 1] = -V1[N - 2];
		Ro1[N - 1] = Ro1[N - 2];

		temp = E1;
		E1 = E;
		E = temp;

		temp = P1;
		P1 = P;
		P = temp;

		temp = V1;
		V1 = V;
		V = temp;

		temp = Ro1;
		Ro1 = Ro;
		Ro = temp;

		t += dt;
		if (loopCounter % 100 == 0)
			cout << loopCounter << " " << t << endl;
	}

}