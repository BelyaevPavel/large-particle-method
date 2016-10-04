#include "LPM.h"


LPM::LPM()
{
	//Gas and particles constant parameters initialisation
	gamma = 1.4;
	R = 287;
	cp = 1005;
	cv = cp / gamma;
	c2 = 710;
	ro01 = 1.21;
	ro02 = 7800;
	lambda1 = 0.026;
	p0 = 0.1e6;
	M = 4.2;
	mu1 = 1.85e-5;
	T0 = 293;
	//Shielding layer constant parameters initialisation
	alpha10 = 1;
	d = 0.001;
	lf = 0.5;
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
	Pr = cp*mu1 / lambda1;
	a = 0.9;
	b = 0.5;
	loopCounter = 0;
	t = 0;
	tstop = 0.01;
	//Grid parameters initialisation
	N = 400;
	Co = 0.01;
	dx = lw / N;
	dt = Co*dx / (sqrt(u1f*u1f + a1*a1));
	//Arrays initialisation
	envVelocity = new double[N];
	U = new double[N];
	UIntermediate = new double[N];
	U1 = new double[N];
	V = new double[N];
	VIntermediate = new double[N];
	V1 = new double[N];
	envE = new double[N];
	gasE = new double[N];
	gasEIntermediate = new double[N];
	gasE1 = new double[N];
	particleE = new double[N];
	particleEIntermediate = new double[N];
	particleE1 = new double[N];
	envPressure = new double[N];
	gasP = new double[N];
	gasPIntermediate = new double[N];
	gasP1 = new double[N];
	particleP = new double[N];
	particlePIntermediate = new double[N];
	particleP1 = new double[N];
	envRo = new double[N];
	gasRo = new double[N];
	gasRo1 = new double[N];
	particleRo = new double[N];
	particleRo1 = new double[N];
	p_count = new double[N];
	p_count1 = new double[N];
	mixtureS = new double[N];
	gasS = new double[N];
	particleS = new double[N];
	//Auxiliary arrays initialisation
	envR = new double[N];
	gasR = new double[N];
	particleR = new double[N];
	envQ = new double[N];
	gasQ = new double[N];
	particleQ = new double[N];
	envF = new double[N];
	gasF = new double[N];
	particleF = new double[N];
	envFi = new double[N];
	gasFi = new double[N];
	particleFi = new double[N];
	alpha = new double[N];
	eta = new double[N];
	Re = new double[N];
	Nu = new double[N];
	Cd = new double[N];
	f = new double[N];
	q = new double[N];
	viscosity1 = new double[N];
	viscosity2 = new double[N];
	gasRightBorderMassFlow = new double[N];
	particleRightBorderMassFlow = new double[N];
	RightBorderParticleFlow = new double[N];

	for (int i = 0; i < N; i++)
	{
		if (i*dx < lf)
		{
			//V[i] = u1f * (i*dx) / lf;
			//xi = pow(1 - (gamma - 1)*(u1f - V[i]) / (2 * a1), 2 / (gamma - 1));
			//envPressure[i] = p1*pow(xi, gamma);
			//gasRo[i] = ro01f*xi;
			envPressure[i] = p1;
			V[i] = u1f;
			gasP[i] = p1;
			gasRo[i] = ro01f;
			U[i] = 0;
			particleRo[i] = 0;
			gasE[i] = envPressure[i] / ((gamma - 1)*gasRo[i]) + V[i] * V[i] / 2.0;
			particleE[i] = 0;
			alpha[i] = 1;
		}
		else
		{
			V[i] = 0;
			U[i] = 0;
			envPressure[i] = 100000;
			gasRo[i] = ro01*alpha10;
			particleRo[i] = ro02 *(1 - alpha10);
			gasE[i] = e0;
			particleE[i] = c2*T0 + U[i] * U[i];
			alpha[i] = alpha10;
		}
		p_count[i] = 6 * (1 - alpha[i]) / (M_PI*d*d*d);
		eta[i] = gasRo[i] * alpha[i] / (gasRo[i] * alpha[i] + particleRo[i] * (1 - alpha[i]));
		gasP[i] = alpha[i] * envPressure[i];
		particleP[i] = (1 - alpha[i]) * envPressure[i];
		envVelocity[i] = eta[i] * V[i] + (1 - eta[i])*U[i];
		envE[i] = eta[i] * (gasP[i] / ((gamma - 1)*gasRo[i]) + 0.5*(envVelocity[i] - V[i])*(envVelocity[i] - V[i])) + (1 - eta[i])* (c2*T0 + 0.5*(envVelocity[i] - U[i])*(envVelocity[i] - U[i]));
		envRo[i] = alpha[i] * gasRo[i] + (1 - alpha[i])*particleRo[i];
	}
}


LPM::~LPM()
{
	delete[] U;
	delete[] UIntermediate;
	delete[] U1;
	delete[] V;
	delete[] VIntermediate;
	delete[] V1;
	delete[] gasE;
	delete[] gasEIntermediate;
	delete[] gasE1;
	delete[] particleE;
	delete[] particleEIntermediate;
	delete[] particleE1;
	delete[] gasP;
	delete[] gasPIntermediate;
	delete[] gasP1;
	delete[] gasRo;
	delete[] gasRo1;
	delete[] particleRo;
	delete[] particleRo1;
	delete[] p_count;
	delete[] p_count1;
	delete[] alpha;
	delete[] eta;
	delete[] Re;
	delete[] Nu;
	delete[] Cd;
	delete[] f;
	delete[] q;
	delete[] viscosity1;
	delete[] viscosity2;
	delete[] gasRightBorderMassFlow;
	delete[] particleRightBorderMassFlow;
	delete[] RightBorderParticleFlow;
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

void LPM::FileOutput(int N, double *Value, double dx, string ValueName, string Aux)
{
	ofstream out(".\\LPM_Output\\"+ValueName + Aux + ".txt");
	for (int i = 0; i < N; i++)
	{
		out <<i*dx<<" "<< Value[i] << endl;
	}
}

void LPM::MainProc()
{
	//Local auxiliary arrays declaration
	double *gasAlpha_05 = new double[N];
	double *particleAlpha_05 = new double[N];
	double *p_count_05 = new double[N];
	double *q_05 = new double[N];
	double *f_05 = new double[N];
	double *gasF_05 = new double[N];
	double *particleF_05 = new double[N];
	double *gasQ_05 = new double[N];
	double *particleQ_05 = new double[N];
	double *gasP_05 = new double[N];
	double *particleP_05 = new double[N];
	double *V_05 = new double[N];
	double *U_05 = new double[N];
	double *dP1 = new double[N];
	double *dP2 = new double[N];
	double *df1 = new double[N];
	double *df2 = new double[N];
	double *dV = new double[N];

	double* temp;
	int n = this->N;
	while (t <= tstop)
	{
		loopCounter++;
		// Artificial viscosity calculation.
		for (int i = 1; i < N - 1; i++)
		{
			if (V[i + 1]<V[i])
				viscosity1[i] = ((gasRo[i + 1] + gasRo[i]) / 2)*(a*sqrt(1.4*(gasP[i + 1] + gasP[i]) / (gasRo[i + 1] + gasRo[i])) + b*(V[i] - V[i + 1]))*(V[i] - V[i + 1]);
			else
				viscosity1[i] = 0;

			if (V[i - 1]>V[i])
				viscosity2[i] = -((gasRo[i - 1] + gasRo[i]) / 2)*(a*sqrt(1.4*(gasP[i - 1] + gasP[i]) / (gasRo[i - 1] + gasRo[i])) - b*(V[i] - V[i - 1]))*(V[i] - V[i - 1]);
			else
				viscosity2[i] = 0;
		}
		viscosity1[0] = viscosity1[1];
		viscosity1[N - 1] = viscosity1[N - 2];
		viscosity2[0] = viscosity2[1];
		viscosity2[N - 1] = viscosity2[N - 2];
		// Auxiliary
		for (int i = 0; i < N; i++)
		{
			Re[i] = gasRo[i] / alpha[i] * fabs(V[i] - U[i])*d / mu1;
			Nu[i] = 2 + 0.6*pow(Pr, 1.0 / 3)*pow(Re[i], 0.5);
			Cd[i] = 0.42*fabs(V[i] - U[i]) + (1 - alpha[i])*mu1 / (gasRo[i] * d)*(24 + 4.4*mu1*sqrt(fabs(V[i] - U[i]))*sqrt(gasRo[i] * d / (alpha[i] * mu1)));
			f[i] = M_PI*d*d*Cd[i] * (gasRo[i] / alpha[i]);
			q[i] = M_PI*d*lambda1*Nu[i] * ((gasE[i] - V[i] * V[i] / 2.0) / cv - (particleE[i] - V[i] * V[i] / 2.0) / c2);
			gasR[i] = (1 - alpha[i])*f[i];
			particleR[i] = -alpha[i] * f[i];
			gasFi[i] = (1 - alpha[i])*q[i];
			particleFi[i] = -alpha[i] * q[i];
			gasF[i] = -0.5*gasRo[i] * (envVelocity[i] - V[i])*(envVelocity[i] - V[i]);
			particleF[i] = -0.5*particleRo[i] * (envVelocity[i] - U[i])*(envVelocity[i] - U[i]);
			envF[i] = -(alpha[i] * gasF[i] + (1 - alpha[i])*particleF[i]);
			gasQ[i] = 0.5*(envVelocity[i] - V[i])*(gasP[i] + gasRo[i] * gasE[i]);
			particleQ[i] = 0.5*(envVelocity[i] - U[i])*(particleP[i] + particleRo[i] * particleE[i]);
			envQ[i] = -(alpha[i] * gasQ[i] + (1 - alpha[i])*particleQ[i]);
		}
		// Intercellar values calc
		for (int i = 0; i < N - 1; i++)
		{
			gasAlpha_05[i] = (alpha[i] + alpha[i + 1]) / 2.0;
			particleAlpha_05[i] = ((1 - alpha[i]) + (1 - alpha[i + 1])) / 2.0;
			p_count_05[i] = (p_count[i] + p_count[i + 1]) / 2.0;
			q_05[i] = (q[i] + q[i + 1]) / 2.0;
			f_05[i] = (f[i] + f[i + 1]) / 2.0;
			gasF_05[i] = (gasF[i] + gasF[i + 1]) / 2.0;
			particleF_05[i] = (particleF[i] + particleF[i + 1]) / 2.0;
			gasQ_05[i] = (gasQ[i] + gasQ[i + 1]) / 2.0;
			particleQ_05[i] = (particleQ[i] + particleQ[i + 1]) / 2.0;
			V_05[i] = (V[i] + V[i + 1]) / 2.0;
			U_05[i] = (U[i] + U[i + 1]) / 2.0;
			gasP_05[i] = (gasP[i] + gasP[i + 1]) / 2.0;
			particleP_05[i] = (particleP[i] + particleP[i + 1]) / 2.0;
		}
		// Intermediate pressure calculation.
		for (int i = 0; i < N - 1; i++)
		{
			gasPIntermediate[i] = (gasP[i + 1] + gasP[i]) / 2.0;
			if (gasPIntermediate[i] < 0)
				cout << loopCounter << "gas Fail." << endl;
			particlePIntermediate[i] = (particleP[i + 1] + particleP[i]) / 2.0;
			if (particlePIntermediate[i] < 0)
				cout << loopCounter << "particle Fail." << endl;
		}
		gasPIntermediate[N - 1] = gasPIntermediate[N - 2];
		particlePIntermediate[N - 1] = particlePIntermediate[N - 2];
		// Pressure and force delta
		//for (int i = 1; i < N - 1; i++)
		//{
		//	dP1[i] = (gasPIntermediate[i] - gasPIntermediate[i - 1] + viscosity1[i] - viscosity2[i])*dt / (dx*(gasRo[i] / alpha[i]));
		//	df1[i] = p_count[i] / gasRo[i] * M_PI*d*d*(gasRo[i] / alpha[i])*Cd[i] * dt / 8.0;

		//	if (p_count[i] <= 1e-50)
		//	{
		//		dP2[i] = 0;
		//		df2[i] = 0;
		//	}
		//	else
		//	{
		//		dP2[i] = 1.0 / (particleRo[i] / (1 - alpha[i]))*(gasPIntermediate[i] - gasPIntermediate[i - 1] + viscosity1[i] - viscosity2[i])*dt / dx;
		//		df2[i] = -p_count[i] / particleRo[i] * M_PI*d*d*(gasRo[i] / alpha[i])*Cd[i] * dt / 8.0;
		//	}
		//}

		//dP1[0] = dP1[1];
		//dP1[N - 1] = dP1[N - 2];

		//dP2[0] = dP2[1];
		//dP2[N - 1] = dP2[N - 2];

		//df1[0] = df1[1];
		//df1[N - 1] = df1[N - 2];

		//df2[0] = df2[1];
		//df2[N - 1] = df2[N - 2];
		// Intermediate velocity calculation.
		for (int i = 1; i < N - 1; i++)
		{
			// dV[i] = ((V[i] - U[i]) - (dP1[i] - dP2[i])) / (1 - (-df1[i] + df2[i]));

			// VIntermediate[i] = V[i] - dP1[i] - df1[i] * dV[i];

			VIntermediate[i] = V[i] - dt / gasRo[i] * (((gasAlpha_05[i] - gasAlpha_05[i - 1])*(gasP[i] + gasF[i]) + alpha[i] * (gasP_05[i] - gasP_05[i - 1] + gasF_05[i] - gasF_05[i - 1])) / dx - alpha[i] * gasR[i]);

			if (p_count[i] > 1e-50)
				UIntermediate[i] = U[i] - dt / particleRo[i] * ((particleP_05[i] - particleP_05[i - 1] + particleF_05[i] - particleF_05[i - 1]) / dx - particleR[i]);
			else
				UIntermediate[i] = 0;
		}
		VIntermediate[0] = VIntermediate[1];
		VIntermediate[N - 1] = -VIntermediate[N - 2];
		UIntermediate[0] = UIntermediate[1];
		UIntermediate[N - 1] = UIntermediate[N - 2];

		for (int i = 1; i < N - 1; i++)
		{
			gasEIntermediate[i] = gasE[i] - dt / gasRo[i] * (1 / dx*((gasAlpha_05[i] - gasAlpha_05[i - 1])*(gasQ[i] + gasF[i] * V[i] + gasP[i] * V[i]) + alpha[i] * ((V_05[i] - V_05[i - 1])*(gasF[i] + gasP[i]) + V[i] * (gasF_05[i] - gasF_05[i - 1] + gasP_05[i] - gasP_05[i - 1]) + (gasQ_05[i] - gasQ_05[i - 1]))) + (alpha[i] * gasFi[i] + gasR[i] * (V[i] + U[i]) / 2.0));
			//gasEIntermediate[i] = gasE[i] + ((particleE[i] - particleEIntermediate[i])*particleRo[i] - ((gasPIntermediate[i] + viscosity1[i])*(gasAlpha_05[i] * (VIntermediate[i] + VIntermediate[i + 1]) / 2 + particleAlpha_05[i] * (UIntermediate[i] + UIntermediate[i + 1]) / 2) - (gasPIntermediate[i - 1] + viscosity2[i])*(gasAlpha_05[i] * (VIntermediate[i] + VIntermediate[i - 1]) / 2 + particleAlpha_05[i] * (UIntermediate[i] + UIntermediate[i - 1]) / 2))*(dt / dx)) / gasRo[i];

			if (gasEIntermediate[i] < 0)
				cout << loopCounter << "Fail." << endl;
			if (p_count[i]>1e-50)
				particleEIntermediate[i] = ((particleE[i] - U[i] * U[i] / 2.0) + U[i] - UIntermediate[i] - dt / ((1 - alpha[i]) * particleRo[i])*(1 / dx*((1 - alpha[i]) * (U_05[i] - U_05[i - 1])*(particleF[i] + particleP[i]) + (particleAlpha_05[i] - particleAlpha_05[i - 1])*(particleQ[i] + particleF[i] * U[i] + particleP[i] * U[i]) + U[i] * (1 - alpha[i]) * (particleF_05[i] - particleF_05[i - 1] + particleP_05[i] - particleP_05[i - 1]) + (1 - alpha[i]) * (particleQ_05[i] - particleQ_05[i - 1])) - particleR[i] * (V[i] + U[i]) / 2.0)) + UIntermediate[i] * UIntermediate[i] * 0.5;
			else
				particleEIntermediate[i] = 0;
			if (particleEIntermediate[i] < 0)
				cout << loopCounter << "Fail." << endl;
		}
		gasEIntermediate[0] = gasEIntermediate[1];
		gasEIntermediate[N - 1] = gasEIntermediate[N - 2];
		particleEIntermediate[0] = particleEIntermediate[1];
		particleEIntermediate[N - 1] = particleEIntermediate[N - 2];
		// Calculation of mass flow through right cell border
		for (int i = 0; i < N - 1; i++)
		{
			gasRightBorderMassFlow[i] = DefineRightBorderMassTranslation(VIntermediate[i], VIntermediate[i + 1], gasRo[i], gasRo[i + 1], dt);

			particleRightBorderMassFlow[i] = DefineRightBorderMassTranslation(UIntermediate[i], UIntermediate[i + 1], particleRo[i], particleRo[i + 1], dt);

			RightBorderParticleFlow[i] = DefineRightBorderMassTranslation(UIntermediate[i], UIntermediate[i + 1], p_count[i], p_count[i + 1], dt);
		}
		RightBorderParticleFlow[N - 1] = 0;
		gasRightBorderMassFlow[N - 1] = gasRightBorderMassFlow[N - 2];
		particleRightBorderMassFlow[N - 1] = 0;
		// Next time step density calculation
		for (int i = 1; i < N; i++)
		{
			gasRo1[i] = gasRo[i] + (gasRightBorderMassFlow[i - 1] - gasRightBorderMassFlow[i]) / dx;

			particleRo1[i] = particleRo[i] + (particleRightBorderMassFlow[i - 1] - particleRightBorderMassFlow[i]) / dx;

			p_count1[i] = p_count[i] + (RightBorderParticleFlow[i - 1] - RightBorderParticleFlow[i]) / dx;

			alpha[i] = 1 - p_count1[i] * M_PI*d*d*d / 6.0;
		}
		gasRo1[0] = gasRo1[1];
		particleRo1[0] = particleRo1[1];
		p_count1[0] = p_count1[1];
		alpha[0] = alpha[1];
		// Final calculations cycle
		for (int i = 1; i < N - 1; i++)
		{
			gasLeftCellVelocityFlow = DefineFlowDirection(i - 1, VIntermediate, 1, VIntermediate);
			gasLeftCellEnergyFlow = DefineFlowDirection(i - 1, VIntermediate, 1, gasEIntermediate);
			gasThisCellEnergyFlow = DefineFlowDirection(i, VIntermediate, 1, gasEIntermediate);
			gasThisCellVelocityFlow = DefineFlowDirection(i, VIntermediate, 1, VIntermediate);

			particleLeftCellVelocityFlow = DefineFlowDirection(i - 1, UIntermediate, 1, UIntermediate);
			particleLeftCellEnergyFlow = DefineFlowDirection(i - 1, UIntermediate, 1, particleEIntermediate);
			particleThisCellEnergyFlow = DefineFlowDirection(i, UIntermediate, 1, particleEIntermediate);
			particleThisCellVelocityFlow = DefineFlowDirection(i, UIntermediate, 1, UIntermediate);

			gasE1[i] = gasRo[i] / gasRo1[i] * gasEIntermediate[i] + (1.0 / (gasRo1[i] * dx))*(gasLeftCellEnergyFlow * gasRightBorderMassFlow[i - 1] - gasThisCellEnergyFlow * gasRightBorderMassFlow[i]);
			V1[i] = gasRo[i] / gasRo1[i] * VIntermediate[i] + (1.0 / (gasRo1[i] * dx))*(gasLeftCellVelocityFlow * gasRightBorderMassFlow[i - 1] - gasThisCellVelocityFlow * gasRightBorderMassFlow[i]);

			if (p_count[i] <= 1e-50)
			{
				particleE1[i] = 0;
				U1[i] = 0;
			}
			else
			{
				particleE1[i] = particleRo[i] / particleRo1[i] * particleEIntermediate[i] + (1.0 / (particleRo1[i] * dx))*(particleLeftCellEnergyFlow * particleRightBorderMassFlow[i - 1] - particleThisCellEnergyFlow * particleRightBorderMassFlow[i]);
				U1[i] = particleRo[i] / particleRo1[i] * UIntermediate[i] + (1.0 / (particleRo1[i] * dx))*(particleLeftCellVelocityFlow * particleRightBorderMassFlow[i - 1] - particleThisCellVelocityFlow * particleRightBorderMassFlow[i]);
			}

			gasP1[i] = gasRo1[i] * (gasE1[i] - V1[i] * V1[i] / 2.0)*0.4;
			particleP1[i] = gasP1[i] / alpha[i] * (1 - alpha[i]);

			//Enviremental values recalculate:
			envPressure[i] = alpha[i] * gasP[i] + (1 - alpha[i])*particleP[i];
			envRo[i] = alpha[i] * gasRo[i] + (1 - alpha[i])*particleRo[i];
			envVelocity[i] = eta[i] * V[i] + (1 - eta[i])*U[i];
			envE[i] = eta[i] * (gasE1[i] + 0.5*(envVelocity[i] - V[i])*(envVelocity[i] - V[i])) + (1 - eta[i])* (particleE[i] + 0.5*(envVelocity[i] - U[i])*(envVelocity[i] - U[i]));
		}

		gasE1[0] = gasE1[1];
		gasE1[N - 1] = gasE1[N - 2];

		particleE1[0] = particleE1[1];
		particleE1[N - 1] = particleE1[N - 2];

		gasP1[0] = gasP1[1];
		gasP1[N - 1] = gasP1[N - 2];

		particleP1[0] = particleP1[1];
		particleP1[N - 1] = particleP1[N - 2];

		V1[0] = V1[1];
		V1[N - 1] = -V1[N - 2];

		U1[0] = U1[1];
		U1[N - 1] = U1[N - 2];

		envE[0] = envE[1];
		envE[N - 1] = envE[N - 2];

		envPressure[0] = envPressure[1];
		envPressure[N - 1] = envPressure[N - 2];

		envRo[0] = envRo[1];
		envRo[N - 1] = envRo[N - 2];

		envVelocity[0] = envVelocity[1];
		envVelocity[N - 1] = envVelocity[N - 2];

		gasRo1[N - 1] = gasRo1[N - 2];

		particleRo1[N - 1] = particleRo1[N - 2];

		p_count1[N - 1] = p_count1[N - 2];

		alpha[N - 1] = alpha[N - 2];

		temp = gasE1;
		gasE1 = gasE;
		gasE = temp;

		temp = particleE1;
		particleE1 = particleE;
		particleE = temp;

		temp = gasP1;
		gasP1 = gasP;
		gasP = temp;

		temp = particleP1;
		particleP1 = particleP;
		particleP = temp;

		temp = V1;
		V1 = V;
		V = temp;

		temp = U1;
		U1 = U;
		U = temp;

		temp = gasRo1;
		gasRo1 = gasRo;
		gasRo = temp;

		temp = particleRo1;
		particleRo1 = particleRo;
		particleRo = temp;

		temp = p_count1;
		p_count1 = p_count;
		p_count = temp;

		if (t > 0.0015)
		{
			FileOutput(N, gasP, dx, "gasP", "t=0.0015");
			break;
		}

		t += dt;
		if (loopCounter % 500 == 0)
			cout << loopCounter << " " << t << endl;
	}

}