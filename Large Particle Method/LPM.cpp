#include "LPM.h"


LPM::LPM()
{
	//Gas and particles constant parameters initialisation
	gamma1 = 1.4;
	gamma2 = 1.5;
	R = 287;
	cp = 1005;
	cv = cp / gamma1;
	c2 = 710;
	ro01 = 1.21;
	ro02 = 1.21;//7800;
	lambda1 = 0.026;
	p0 = 0.1e6;
	M = 4.2;
	mu1 = 1.85e-5;
	T0 = 293;
	//Shielding layer constant parameters initialisation
	alpha10 = 0.5;
	d = 0.001;
	lf = 0.5;
	ls = 3;
	lw = 3.5;
	//Gas calculated parameters initialisation
	p1 = (1 + 2 * gamma1*(M*M - 1) / (gamma1 + 1))*p0;//(1 + 2 * gamma1 / (gamma1 + 1)*M*M - 1)*p0;
	ro01f = (gamma1 + 1)*M*M / (2 + (gamma1 - 1)*M*M)*ro01;
	a10 = sqrt(gamma1*p0 / ro01);
	a1 = sqrt(gamma1*p1 / ro01f);
	u1f = 2 / (gamma1 + 1)*(M - 1 / M)*a10;
	e0 = p0 / (ro01*(gamma1 - 1));
	e1 = p1 / (ro01f*(gamma1 - 1)) + 0.5*u1f*u1f;
	//Auxiliary variables initialisation
	Pr = cp*mu1 / lambda1;
	a = 0.9;
	b = 0.5;
	a12 = 0.5;
	a21 = 0.5;
	b12 = 0.009;
	b21 = 0.002;
	d12 = 0.00005;//0.05;
	d21 = 0.00005;//0.05;
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
	sigma1 = new double[N];
	sigma2 = new double[N];
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
	gasAlpha = new double[N];
	gasAlphaIntermediate = new double[N];
	gasAlpha1 = new double[N];
	particleAlpha = new double[N];
	particleAlphaIntermediate = new double[N];
	particleAlpha1 = new double[N];
	gasEta = new double[N];
	particleEta = new double[N];
	Re = new double[N];
	Nu = new double[N];
	Cd = new double[N];
	f = new double[N];
	q = new double[N];
	gas1f = new double[N];
	gas1q = new double[N];
	gas2f = new double[N];
	gas2q = new double[N];
	viscosity1 = new double[N];
	viscosity2 = new double[N];
	gasRightBorderMassFlow = new double[N];
	particleRightBorderMassFlow = new double[N];
	RightBorderGasAlphaFlow = new double[N];

	for (int i = 0; i < N; i++)
	{
		particleAlpha1[i] = 0;
		particleE1[i] = 0;
		particleP1[i] = 0;
		particleRo1[i] = 0;
		gasAlpha1[i] = 0;
		gasE1[i] = 0;
		gasP1[i] = 0;
		gasRo1[i] = 0;
		V1[i] = 0;
		if (i*dx < lf)//(i*dx < lf)
		{
			//V[i] = u1f * (i*dx) / lf;
			//xi = pow(1 - (gamma1 - 1)*(u1f - V[i]) / (2 * a1), 2 / (gamma1 - 1));
			//envPressure[i] = p1*pow(xi, gamma1);
			//gasRo[i] = ro01f*xi;
			V[i] = u1f;
			gasP[i] = p1;
			particleP[i] = p1;// p1;
			gasRo[i] = ro01f;//ro01f*alpha10;
			U[i] = u1f;
			particleRo[i] = ro01f;//ro01f*(1 - alpha10);
			gasE[i] = gasP[i] / ((gamma1 - 1)*gasRo[i]) + V[i] * V[i] / 2.0;
			particleE[i] = particleP[i] / ((gamma2 - 1)*particleRo[i]) + U[i] * U[i] / 2.0;
			gasAlpha[i] = alpha10;
			particleAlpha[i] = (1 - alpha10);
		}
		else
		{
			V[i] = 0;
			U[i] = 0;
			gasP[i] = p0;// 100000;
			particleP[i] = p0;// 100000;
			gasRo[i] = ro01;// *alpha10;
			particleRo[i] = ro02;// *(1 - alpha10);
			gasE[i] = gasP[i] / ((gamma1 - 1)*gasRo[i]) + V[i] * V[i] / 2.0;
			particleE[i] = particleP[i] / ((gamma1 - 1)*particleRo[i]) + U[i] * U[i] / 2.0;//c2*T0 + U[i] * U[i];
			gasAlpha[i] = alpha10;
			particleAlpha[i] = (1 - alpha10);
		}
		p_count[i] = 6 * (1 - gasAlpha[i]) / (M_PI*d*d*d);
		gasEta[i] = gasRo[i] * gasAlpha[i] / (gasRo[i] * gasAlpha[i] + particleRo[i] * particleAlpha[i]);
		particleEta[i] = particleRo[i] * particleAlpha[i] / (gasRo[i] * gasAlpha[i] + particleRo[i] * particleAlpha[i]);
		envPressure[i] = gasAlpha[i] * gasP[i] + particleAlpha[i] * envPressure[i];
		envVelocity[i] = gasEta[i] * V[i] + (1 - gasEta[i])*U[i];
		envE[i] = gasEta[i] * (gasP[i] / ((gamma1 - 1)*gasRo[i]) + 0.5*(envVelocity[i] - V[i])*(envVelocity[i] - V[i])) + (1 - gasEta[i])* (c2*T0 + 0.5*(envVelocity[i] - U[i])*(envVelocity[i] - U[i]));
		envRo[i] = gasAlpha[i] * gasRo[i] + particleAlpha[i] * particleRo[i];
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
	delete[] gasAlpha;
	delete[] gasEta;
	delete[] Re;
	delete[] Nu;
	delete[] Cd;
	delete[] f;
	delete[] q;
	delete[] viscosity1;
	delete[] viscosity2;
	delete[] gasRightBorderMassFlow;
	delete[] particleRightBorderMassFlow;
	delete[] RightBorderGasAlphaFlow;
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
	ofstream out(".\\LPM_Output\\" + ValueName + Aux + ".txt");
	for (int i = 0; i < N; i++)
	{
		out << i*dx << " " << Value[i] << endl;
	}
}

void LPM::MainProc(int MixtureType)
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
			//if (V[i + 1]<V[i])
			//	viscosity1[i] = ((gasRo[i + 1] + gasRo[i]) / 2)*(a*sqrt(1.4*(gasP[i + 1] + gasP[i]) / (gasRo[i + 1] + gasRo[i])) + b*(V[i] - V[i + 1]))*(V[i] - V[i + 1]);
			//else
			viscosity1[i] = 0;

			//if (V[i - 1]>V[i])
			//	viscosity2[i] = -((gasRo[i - 1] + gasRo[i]) / 2)*(a*sqrt(1.4*(gasP[i - 1] + gasP[i]) / (gasRo[i - 1] + gasRo[i])) - b*(V[i] - V[i - 1]))*(V[i] - V[i - 1]);
			//else
			viscosity2[i] = 0;
		}
		viscosity1[0] = viscosity1[1];
		viscosity1[N - 1] = viscosity1[N - 2];
		viscosity2[0] = viscosity2[1];
		viscosity2[N - 1] = viscosity2[N - 2];
		// Auxiliary
		for (int i = 0; i < N; i++)
		{
			switch (MixtureType)
			{
				case 1:
					gas1f[i] = -0.5*gasRo[i] * (envVelocity[i] - V[i])*(envVelocity[i] - V[i]);
					gas2f[i] = -0.5*particleRo[i] * (envVelocity[i] - U[i])*(envVelocity[i] - U[i]);
					gas1q[i] = 0.5*(envVelocity[i] - V[i])*(gasP[i] + gasRo[i] * gasE[i]);
					gas2q[i] = 0.5*(envVelocity[i] - U[i])*(particleP[i] + particleRo[i] * particleE[i]);
					gasR[i] = particleAlpha[i] * a21*(U[i] - V[i]);
					particleR[i] = gasAlpha[i] * a12*(V[i] - U[i]);
					gasFi[i] = (1 - gasAlpha[i])*b21*(particleP[i] - gasP[i]);
					particleFi[i] = gasAlpha[i] * b12*(gasP[i] - particleP[i]);
					gasF[i] = gas1f[i];
					particleF[i] = gas2f[i];
					envF[i] = -(gasAlpha[i] * gasF[i] + (1 - gasAlpha[i])*particleF[i]);
					gasQ[i] = gas1q[i];
					particleQ[i] = gas2q[i];
					envQ[i] = -(gasAlpha[i] * gasQ[i] + (1 - gasAlpha[i])*particleQ[i]);
					sigma1[i] = gasAlpha[i] > 1e-50 ? particleAlpha[i] * d12*b12*(gasP[i] - particleP[i]) : 0;
					sigma2[i] = particleAlpha[i] > 1e-50 ? gasAlpha[i] * d21*b21*(particleP[i] - gasP[i]) : 0;
					break;
				case 2:
					Re[i] = gasRo[i] / gasAlpha[i] * fabs(V[i] - U[i])*d / mu1;
					Nu[i] = 2 + 0.6*pow(Pr, 1.0 / 3)*pow(Re[i], 0.5);
					Cd[i] = 0.42*fabs(V[i] - U[i]) + (1 - gasAlpha[i])*mu1 / (gasRo[i] * d)*(24 + 4.4*mu1*sqrt(fabs(V[i] - U[i]))*sqrt(gasRo[i] * d / (gasAlpha[i] * mu1)));
					f[i] = M_PI*d*d*Cd[i] * (gasRo[i] / gasAlpha[i]);
					q[i] = M_PI*d*lambda1*Nu[i] * ((gasE[i] - V[i] * V[i] / 2.0) / cv - (particleE[i] - V[i] * V[i] / 2.0) / c2);
					gasR[i] = (1 - gasAlpha[i])*f[i];
					particleR[i] = -gasAlpha[i] * f[i];
					gasFi[i] = (1 - gasAlpha[i])*q[i];
					particleFi[i] = -gasAlpha[i] * q[i];
					gasF[i] = -0.5*gasRo[i] * (envVelocity[i] - V[i])*(envVelocity[i] - V[i]);
					particleF[i] = -0.5*particleRo[i] * (envVelocity[i] - U[i])*(envVelocity[i] - U[i]);
					envF[i] = -(gasAlpha[i] * gasF[i] + (1 - gasAlpha[i])*particleF[i]);
					gasQ[i] = 0.5*(envVelocity[i] - V[i])*(gasP[i] + gasRo[i] * gasE[i]);
					particleQ[i] = 0.5*(envVelocity[i] - U[i])*(particleP[i] + particleRo[i] * particleE[i]);
					envQ[i] = -(gasAlpha[i] * gasQ[i] + (1 - gasAlpha[i])*particleQ[i]);
					break;
				default:
					break;
			}
		}
		// Intercellar values calc
		for (int i = 0; i < N - 1; i++)
		{
			gasAlpha_05[i] = (gasAlpha[i] + gasAlpha[i + 1]) / 2.0;
			particleAlpha_05[i] = (particleAlpha[i] + particleAlpha[i + 1]) / 2.0;
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
		gasAlpha_05[N - 1] = gasAlpha_05[N - 2];
		particleAlpha_05[N - 1] = particleAlpha_05[N - 2];
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
		//	dP1[i] = (gasPIntermediate[i] - gasPIntermediate[i - 1] + viscosity1[i] - viscosity2[i])*dt / (dx*(gasRo[i] / gasAlpha[i]));
		//	df1[i] = p_count[i] / gasRo[i] * M_PI*d*d*(gasRo[i] / gasAlpha[i])*Cd[i] * dt / 8.0;

		//	if (p_count[i] <= 1e-50)
		//	{
		//		dP2[i] = 0;
		//		df2[i] = 0;
		//	}
		//	else
		//	{
		//		dP2[i] = 1.0 / (particleRo[i] / (1 - gasAlpha[i]))*(gasPIntermediate[i] - gasPIntermediate[i - 1] + viscosity1[i] - viscosity2[i])*dt / dx;
		//		df2[i] = -p_count[i] / particleRo[i] * M_PI*d*d*(gasRo[i] / gasAlpha[i])*Cd[i] * dt / 8.0;
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

			//gasAlphaIntermediate[i] = gasAlpha[i] - dt*(envVelocity[i] * (gasAlpha[i+1] - gasAlpha[i - 1]) / (2*dx) + gasAlpha[i] * sigma1[i]);

			//particleAlphaIntermediate[i] = particleAlpha[i] - dt*(envVelocity[i] * (particleAlpha[i+1] - particleAlpha[i - 1]) / (2*dx) + particleAlpha[i] * sigma2[i]);

			if (gasAlpha[i] > 1e-50)
				VIntermediate[i] = V[i] - dt / gasRo[i] * ((gasP_05[i] - gasP_05[i - 1] + gasF_05[i] - gasF_05[i - 1]) / dx - gasR[i]);//V[i] - dt / gasRo[i] * (((gasAlpha_05[i] - gasAlpha_05[i - 1])*(gasP[i] + gasF[i]) + gasAlpha[i] * (gasP_05[i] - gasP_05[i - 1] + gasF_05[i] - gasF_05[i - 1])) / dx - gasAlpha[i] * gasR[i]);
			else
				VIntermediate[i] = 0;

			if (particleAlpha[i] > 1e-50)
				UIntermediate[i] = U[i] - dt / particleRo[i] * ((particleP_05[i] - particleP_05[i - 1] + particleF_05[i] - particleF_05[i - 1]) / dx - particleR[i]);
			else
				UIntermediate[i] = 0;
		}
		//gasAlphaIntermediate[0] = gasAlphaIntermediate[1];
		//gasAlphaIntermediate[N - 1] = gasAlphaIntermediate[N - 2];
		//particleAlphaIntermediate[0] = particleAlphaIntermediate[1];
		//particleAlphaIntermediate[N - 1] = particleAlphaIntermediate[N - 2];
		VIntermediate[0] = VIntermediate[1];
		VIntermediate[N - 1] = -VIntermediate[N - 2];
		UIntermediate[0] = UIntermediate[1];
		UIntermediate[N - 1] = -UIntermediate[N - 2];

		for (int i = 1; i < N - 1; i++)
		{
			if (gasAlpha[i] > 1e-50)
			{
				gasEIntermediate[i] = gasE[i] - dt / gasRo[i] * (1 / dx*((V_05[i] - V_05[i - 1])*(gasF[i] + gasP[i]) + V[i] * (gasF_05[i] - gasF_05[i - 1] + gasP_05[i] - gasP_05[i - 1]) + (gasQ_05[i] - gasQ_05[i - 1])) + (gasFi[i] + gasAlpha[i] * gasR[i] * (V[i] + U[i]) / 2.0));
				if (gasEIntermediate[i] < 0)
					cout << loopCounter << "Fail." << endl;
			}
			else
				gasEIntermediate[i] = 0;
			//gasEIntermediate[i] = gasE[i] + ((particleE[i] - particleEIntermediate[i])*particleRo[i] - ((gasPIntermediate[i] + viscosity1[i])*(gasAlpha_05[i] * (VIntermediate[i] + VIntermediate[i + 1]) / 2 + particleAlpha_05[i] * (UIntermediate[i] + UIntermediate[i + 1]) / 2) - (gasPIntermediate[i - 1] + viscosity2[i])*(gasAlpha_05[i] * (VIntermediate[i] + VIntermediate[i - 1]) / 2 + particleAlpha_05[i] * (UIntermediate[i] + UIntermediate[i - 1]) / 2))*(dt / dx)) / gasRo[i];

			if (particleAlpha[i] > 1e-50)
			{
				particleEIntermediate[i] = particleE[i] - dt / particleRo[i] * (1 / dx*((U_05[i] - U_05[i - 1])*(particleF[i] + particleP[i]) + U[i] * (particleF_05[i] - particleF_05[i - 1] + particleP_05[i] - particleP_05[i - 1]) + (particleQ_05[i] - particleQ_05[i - 1])) + (particleFi[i] + particleAlpha[i] * particleR[i] * (V[i] + U[i]) / 2.0));
				if (particleEIntermediate[i] < 0)
					cout << loopCounter << "Fail." << endl;
			}
			else
				particleEIntermediate[i] = 0;
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

			RightBorderGasAlphaFlow[i] = DefineRightBorderMassTranslation(VIntermediate[i], VIntermediate[i + 1], gasAlpha[i], gasAlpha[i + 1], dt);
		}
		RightBorderGasAlphaFlow[N - 1] = 0;
		gasRightBorderMassFlow[N - 1] = 0;// gasRightBorderMassFlow[N - 2];
		particleRightBorderMassFlow[N - 1] = 0;
		// Next time step density calculation
		for (int i = 1; i < N - 1; i++)
		{
			if (gasAlpha[i] > 1e-50)
				gasRo1[i] = gasRo[i] + (gasRightBorderMassFlow[i - 1] - gasRightBorderMassFlow[i]) / dx;
			else
				gasRo1[i] = 0;

			if (particleAlpha[i] > 1e-50)
				particleRo1[i] = particleRo[i] + (particleRightBorderMassFlow[i - 1] - particleRightBorderMassFlow[i]) / dx;
			else
				particleRo1[i] = 0;

			if (gasAlpha[i] > 1e-50)
				gasAlpha1[i] = gasAlpha[i] + (RightBorderGasAlphaFlow[i - 1] - RightBorderGasAlphaFlow[i]) / dx - sigma1[i];
			else
				gasAlpha1[i] = 0;

			if (particleAlpha[i] > 1e-50)
				particleAlpha1[i] = particleAlpha[i] + (DefineRightBorderMassTranslation(UIntermediate[i - 1], UIntermediate[i], particleAlpha[i - 1], particleAlpha[i], dt) - DefineRightBorderMassTranslation(UIntermediate[i], UIntermediate[i + 1], particleAlpha[i], particleAlpha[i + 1], dt)) / dx - sigma2[i];
			else
				particleAlpha1[i] = 0;

			//gasAlpha[i] = 1 - p_count1[i] * M_PI*d*d*d / 6.0;
			gasEta[i] = gasRo1[i] * gasAlpha1[i] / (gasRo1[i] * gasAlpha1[i] + particleRo1[i] * particleAlpha1[i]);
			particleEta[i] = particleRo1[i] * particleAlpha1[i] / (gasRo1[i] * gasAlpha1[i] + particleRo1[i] * particleAlpha1[i]);

		}
		gasEta[0] = gasEta[1];
		gasEta[N - 1] = gasEta[N - 2];
		particleEta[0] = particleEta[1];
		particleEta[N - 1] = particleEta[N - 2];
		gasRo1[0] = gasRo1[1];
		particleRo1[0] = particleRo1[1];
		gasRo1[N - 1] = gasRo1[N - 2];
		particleRo1[N - 1] = particleRo1[N - 2];
		p_count1[0] = p_count1[1];
		gasAlpha1[0] = gasAlpha1[1];
		gasAlpha1[N - 1] = gasAlpha1[N - 2];
		particleAlpha1[0] = particleAlpha1[1];
		particleAlpha1[N - 1] = particleAlpha1[N - 2];
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

			if (gasAlpha[i] <= 1e-50)
			{
				gasE1[i] = 0;
				V1[i] = 0;
				gasP1[i] = 0;
			}
			else
			{
				gasE1[i] = gasRo[i] / gasRo1[i] * gasEIntermediate[i] + (1.0 / (gasRo1[i] * dx))*(gasLeftCellEnergyFlow * gasRightBorderMassFlow[i - 1] - gasThisCellEnergyFlow * gasRightBorderMassFlow[i]);
				V1[i] = gasRo[i] / gasRo1[i] * VIntermediate[i] + (1.0 / (gasRo1[i] * dx))*(gasLeftCellVelocityFlow * gasRightBorderMassFlow[i - 1] - gasThisCellVelocityFlow * gasRightBorderMassFlow[i]);
				gasP1[i] = gasRo1[i] * (gasE1[i] - V1[i] * V1[i] / 2.0)*0.4;
			}
			//gasAlpha1[i] = gasRo[i] / gasRo1[i] * gasAlphaIntermediate[i] + (1.0 / (gasRo1[i] * dx))*(DefineFlowDirection(i - 1, VIntermediate, 1, gasAlphaIntermediate) * gasRightBorderMassFlow[i - 1] - DefineFlowDirection(i, VIntermediate, 1, gasAlphaIntermediate) * gasRightBorderMassFlow[i]);

			if (particleAlpha[i] <= 1e-50)
			{
				particleE1[i] = 0;
				U1[i] = 0;
				particleP1[i] = 0;
			}
			else
			{
				particleE1[i] = particleRo[i] / particleRo1[i] * particleEIntermediate[i] + (1.0 / (particleRo1[i] * dx))*(particleLeftCellEnergyFlow * particleRightBorderMassFlow[i - 1] - particleThisCellEnergyFlow * particleRightBorderMassFlow[i]);
				U1[i] = particleRo[i] / particleRo1[i] * UIntermediate[i] + (1.0 / (particleRo1[i] * dx))*(particleLeftCellVelocityFlow * particleRightBorderMassFlow[i - 1] - particleThisCellVelocityFlow * particleRightBorderMassFlow[i]);
				particleP1[i] = particleRo1[i] * (particleE1[i] - U1[i] * U1[i] / 2.0)*0.4;
			}


			//Enviremental values recalculate:
			envPressure[i] = gasAlpha1[i] * gasP1[i] + particleAlpha1[i] * particleP1[i];
			envRo[i] = gasAlpha1[i] * gasRo[i] + particleAlpha1[i] * particleRo1[i];
			envVelocity[i] = gasEta[i] * V1[i] + particleEta[i] * U1[i];
			envE[i] = gasEta[i] * (gasE1[i] + 0.5*(envVelocity[i] - V1[i])*(envVelocity[i] - V[i])) + particleEta[i] * (particleE1[i] + 0.5*(envVelocity[i] - U1[i])*(envVelocity[i] - U1[i]));
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
		U1[N - 1] = -U1[N - 2];

		envE[0] = envE[1];
		envE[N - 1] = envE[N - 2];

		envPressure[0] = envPressure[1];
		envPressure[N - 1] = envPressure[N - 2];

		envRo[0] = envRo[1];
		envRo[N - 1] = envRo[N - 2];

		envVelocity[0] = envVelocity[1];
		envVelocity[N - 1] = gasEta[N - 1] * V1[N - 1] + (1 - gasEta[N - 1])*U1[N - 1];;

		gasRo1[N - 1] = gasRo1[N - 2];

		particleRo1[N - 1] = particleRo1[N - 2];

		p_count1[N - 1] = p_count1[N - 2];

		//gasAlpha1[0] = gasAlpha1[1];
		//gasAlpha1[N - 1] = gasAlpha1[N - 2];


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

		temp = gasAlpha1;
		gasAlpha1 = gasAlpha;
		gasAlpha = temp;

		temp = particleAlpha1;
		particleAlpha1 = particleAlpha;
		particleAlpha = temp;
		//if (t > 0.0015)
		//{
		//	FileOutput(N, gasP, dx, "gasP", "t=0.0015");
		//	break;
		//}

		t += dt;
		if (loopCounter % 1000 == 0)
			cout << loopCounter << " " << t << endl;
	}

}