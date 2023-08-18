#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"
#include <complex.h>
#include <fftw3.h>
#include "stuff.h"

void create_plan(void);
void destroy_plan(void);
void Right(fftw_complex *rhsR, fftw_complex *rhsV);
void Shift(fftw_complex *rhsR, fftw_complex *rhsV, double coeffRK);
void write_evolution(void);

void calc_steepness(void){
	int i;
	
	
	dYk[0] = 0.0 + I * 0.0;
	
	for(i = 1; i < N/2; i++){
		dYk[i] = I * ( ( (double)i ) * Dk ) * Yk[i];
		dYk[N-i] = conj( dYk[i] );		
	}
	
	dYk[0] = 0.0 + I * 0.0;
	dYk[N/2] = 0.0 + I * 0.0;

	fftw_execute(dYk_2_dYu_);

	
	dEtaMax = 0.0;
	dEtaMin = 0.0;
	for(i = 0; i < N; i++){
		Steepness_u[i] = creal(dYu[i]) / (1.0 + creal(dXTu[i]) );

		if( dEtaMax < creal(Steepness_u[i]) ){
			dEtaMax = creal(Steepness_u[i]);
		}

		if( dEtaMin > creal(Steepness_u[i]) ){
			dEtaMin = creal(Steepness_u[i]);
		}
	}

	// fftw_execute(Steepness_u_2_Steepness_k);

	// // #pragma omp parallel for private(i) shared(Steepness_k1)
	// for(i = 0; i < N; i++){
	// 	Steepness_k[i] = Steepness_k[i] * norm;
	// }

}


void write_restart(void){
	int number_of_numbers;
	char file_name[40];
	sprintf(file_name, "./restart_R/%d_restart", nstep_Restart);
	if( (restart = fopen (file_name,"w+b"))==NULL )	{
		printf("Cannot open the file restart \n");
	}
	else{

	}
	number_of_numbers = fwrite(&RkCur[0], sizeof(fftw_complex), N, restart);
	fclose(restart);

	sprintf(file_name, "./restart_V/%d_restart", nstep_Restart);
	if( (restart = fopen (file_name,"w+b"))==NULL )	{
		printf("Cannot open the file restart \n");
	}
	else{

	}
	number_of_numbers = fwrite(&VkCur[0], sizeof(fftw_complex), N, restart);
	fclose(restart);
}
double NLSE_soliton(double x){
	return lNLS/sqrt(2.0)/pow((double)k0*Dk, 2.0)/cosh(lNLS*(x- Length/2.0 ))*cos( ((double)k0)*Dk*(x-Length/2.0) );
}
double FindEtax(double x){
	int k;
	fftw_complex etaxCur;
	// etaxCur = 0.0 + I*0.0;
	// for(k=0; k<N; k++){
	// 	etaxCur = etaxCur + Etak[k]*cexp(I*Dk*((double)k)*x);
	// }
	etaxCur = Etak[0];
	for(k=1; k<N/2; k++){
		etaxCur = etaxCur + Etak[k]*cexp(I*Dk*((double)k)*x) + conj(Etak[k]*cexp(I*Dk*((double)k)*x));
		// etaxCur = etaxCur + Etak[N-k]*cexp(-I*Dk*((double)(N-k))*x);
	}	
	// printf("%e   %e   %e\n", x, creal(etaxCur), cimag(etaxCur));
	return creal(etaxCur);
}

void find_z(void){
	// double minErr = 3.5e-2;
	// double minErr = 1.e-16;
	double minErr = 1e-15;
	double error = 1.e6;
	int u, i;
	double uCur;
	double Tek;


	for(u = 0; u < N; u++){
		XTu[u] = 0.0;
		Yu[u] = 0.0;
		// uCur = du*((double)u);
		// Tek = FindEtax(uCur + XTu[u]);
	}
	while(error > minErr){
		error = 0.0;
		for(u = 0; u < N; u++){
			uCur = du * ( (double)u );
			// YuCur[u] = NLSE_soliton(uCur + XTu[u]) + I*0.0;
			YuCur[u] = FindEtax(uCur + creal(XTu[u])) + I*0.0;
			error += fabs(creal(YuCur[u]) - creal(Yu[u]) );
			Yu[u] = YuCur[u];
			// printf("%d   %e\n", u, error);
		}
		error = error*norm;
		

		fftw_execute(Yu_2_Yk);
		XTk[0] = 0.0 + I*0.0;
		Yk[0] = Yk[0]*norm;
		for(i=1; i<N/2; i++){
			Yk[i] = Yk[i]*norm;
			Yk[N-i] = Yk[N-i]*norm;
			XTk[i] = -I*Yk[i];
			XTk[N-i] = I*Yk[N-i];
			// printf("%d   %e   %e\n", i, creal(Yk[i]), cimag(Yk[i]) );
		}
		XTk[N/2] = 0.0 + I*0.0;
		Yk[N/2] = 0.0 + I*0.0;
		fftw_execute(XTk_2_XTu_);
		printf("error = %e   MinError = %e\n", error, minErr);
	
	}

/*	FILE *pfile;
	pfile = fopen("y_k.dat", "w");

	for (int j = 0; j < N; j++)
	{
		fprintf(pfile, "%d\t%e\t%e\t%e\n", j, cabs(Yk[j]), cabs(XTk[j]), cabs(Etak[j]));
	}

	fclose(pfile);*/

	//fftw_execute(Etak_2_Etax);


/*	FILE *eta;
	eta = fopen("nlse_soliton.dat", "w");
	for(u=0; u<N; u++){
		uCur = du*((double)u);
		fprintf(eta, "%e   %e   %e   %e 	%e\n", uCur, creal(Etax[u]), uCur + creal(XTu[u]), creal(Yu[u]), creal(XTu[u]));
	}
	fclose(eta);*/

	for(i=1; i<N/2; i++){
		dZk[i] = 0.0;//
		dZk[N-i] = -I*((double)i)*Dk*(XTk[N-i] + I*Yk[N-i]);
		// dZk[N-i] = I*Diff[N-i]*(XTk[N-i] + I*Yk[N-i]);
	}

	dZk[0] = 1.0 + I*0.0;
	dZk[N/2] = 0.0 + I*0.0;
}

// void find_z(void){
// 	double minErr = 1.7e-19;
// 	double error = 1.e6;
// 	int u, i;
// 	double uCur;
// 	for(u=0; u<N; u++){
// 		XTu[u] = 0.0;
// 		Yu[u] = 0.0;
// 	}
// 	while(error>minErr){
// 		error = 0.0;
// 		for(u=0; u<N; u++){
// 			uCur = du*((double)u);
// 			YuCur[u] = NLSE_soliton(uCur + XTu[u]) + I*0.0;
// 			error += fabs(creal(YuCur[u]) - creal(Yu[u]) );
// 			Yu[u] = YuCur[u];
// 			// printf("%d   %e\n", u, error);
// 		}
// 		error = error*norm;
		

// 		fftw_execute(Yu_2_Yk);
// 		XTk[0] = 0.0 + I*0.0;
// 		Yk[0] = Yk[0]*norm;
// 		for(i=1; i<N/2; i++){
// 			Yk[i] = Yk[i]*norm;
// 			Yk[N-i] = Yk[N-i]*norm;
// 			XTk[i] = -I*Yk[i];
// 			XTk[N-i] = I*Yk[N-i];
// 			// printf("%d   %e   %e\n", i, creal(Yk[i]), cimag(Yk[i]) );
// 		}
// 		XTk[N/2] = 0.0 + I*0.0;
// 		Yk[N/2] = 0.0 + I*0.0;
// 		fftw_execute(XTk_2_XTu_);
// 		printf("error = %e   MinError = %e\n", error, minErr);
	
//	}
/*	FILE *eta;
	eta = fopen("nlse_soliton.dat", "w");
	for(u=0; u<N; u++){
		uCur = du*((double)u);
		fprintf(eta, "%e   %e   %e   %e\n", uCur, NLSE_soliton(uCur), uCur + creal(XTu[u]), creal(Yu[u]));
	}
	fclose(eta);*/
/*	for(i=1; i<N/2; i++){
		dZk[i] = 0.0;//
		dZk[N-i] = -I*((double)i)*Dk*(XTk[N-i] + I*Yk[N-i]);
		// dZk[N-i] = I*Diff[N-i]*(XTk[N-i] + I*Yk[N-i]);
	}
	dZk[0] = 1.0 + I*0.0;
	dZk[N/2] = 0.0 + I*0.0;
}*/

void find_R(void){
	int i, j;
	fftw_complex sum;

	RkCur[0] = 1.0;
	
	for(i=1; i<N/2; i++){
		sum = 0.0 + I*0.0;
		for(j=1; j<i; j++){
			sum = sum + dZk[N-j]*RkCur[N-(i-j)];
		}
		RkCur[N-i] = - dZk[N-i] - sum;
	}	
	for(i=1; i<=N/2; i++){
		RkCur[i] = 0.0 + I*0.0;
	}

}
void read_RV(void){
	int number_of_numbers;
	char file_name[40];
	sprintf(file_name, "./restart_R/%d_restart", contin);
	if( (data = fopen (file_name,"rb"))==NULL )
	{
		printf("Cannot open the file restart_R  %d \n", contin);
	}
	else{
		printf("Open the file restart_R  %d \n", contin);
	}
	number_of_numbers = fread(&RkCur[0], sizeof(fftw_complex), N, data);
	fclose(data);

	sprintf(file_name, "./restart_V/%d_restart", contin);
	if( (data = fopen (file_name,"rb"))==NULL )
	{
		printf("Cannot open the file restart_V  %d \n", contin);
	}
	else{
		printf("Open the file restart_V  %d \n", contin);
	}
	number_of_numbers = fread(&VkCur[0], sizeof(fftw_complex), N, data);
	fclose(data);
}
void read_etak(void){
	int number_of_numbers;
	char file_name[40];
	sprintf(file_name, "./data_etafin");
	if( (data = fopen (file_name,"rb"))==NULL )
	{
		printf("Cannot open the file data_etafin\n");
	}
	else{
		printf("Open the file data_etafin\n");
	}
	number_of_numbers = fread(&Etak[0], sizeof(fftw_complex), N, data);
	fclose(data);

}
void read_psik(void){
	int number_of_numbers;
	char file_name[40];
	sprintf(file_name, "./data_psifin");
	if( (data = fopen (file_name,"rb"))==NULL )
	{
		printf("Cannot open the file data_psifin\n");
	}
	else{
		printf("Open the file data_psifin\n");
	}
	number_of_numbers = fread(&PsikIni[0], sizeof(fftw_complex), N, data);
	fclose(data);

}

void find_V(){
	int u;
	int k;
	double xCur;

	FILE *pfile;

// find Psi(u)	
	fftw_complex PsiuIniSum;
	for(u=0; u<N; u++){
		// uCur = du*((double)u);
		xCur = du*((double)u) + XTu[u];
		PsiuIniSum = PsikIni[0];
		for(k=1; k<N/2; k++){
			PsiuIniSum = PsiuIniSum + PsikIni[k]*cexp(I*Dk*((double)k)*xCur) + conj( PsikIni[k]*cexp(I*Dk*((double)k)*xCur) );
		}
		Psiu[u] = creal(PsiuIniSum);
	}
// find Psi(u)
	fftw_execute(Psiu_2_Psik);

	for (int j = 0; j < N; j++)
	{
		Psik[j] = Psik[j]/N;
	}

	for(k=0; k<=N/2; k++){
		dPhik[k] = 0.0 + I*0.0;
	}
	for(k=1; k<N/2; k++){
		dPhik[N-k] = -2.0*I*Dk*((double)k)*Psik[N-k];
	}
	fftw_execute(dPhik_2_dPhiu_);
	fftw_execute(RkCur_2_RuCur_);

	for(u=0; u<N; u++){
		VuCur[u] = I*RuCur[u]*dPhiu[u];
	}

	fftw_execute(VuCur_2_VkCur);

	for(k=0; k<=N/2; k++){
		VkCur[k] = 0.0 + I*0.0;
	}
	for(k=1; k<N/2; k++){
		VkCur[N-k] = VkCur[N-k]*norm;
	}
}

void initial(void){
	int k;
	double tauMax;
	double phase;
	int rnd_num;

	EtaMax = 0.0;
	EtaMin = 0.0;
	Length = 10000.0;
	// Length = 2.0 * M_PI;
	
	// Tp = 1273.239545;
	// x0 = 5250.0;
	Period_calc_pos = Tp / ( tau * freq_ncalc_pos );
	printf("Period_calc_pos = %d\n", Period_calc_pos);

	du = Length/((double)N);
	Dk = 2.0 * M_PI/Length;

	tauMax = du * sqrt(Dk/grav);

	printf("tauMax = %e   tau = %e\n", tauMax, tau);

	// Vgr = 0.5 * sqrt( grav / (100.0 * Dk) );

	// Vgr = 6.466103;

	printf("FinTime = %f\n", FinTime);



	norm = 1.0/((double)N);

	for(k = 0; k < N/2; k++){
		Diff[k]	= Dk * ( (double)k );
	}
	Diff[N/2] = 0.0 + I * 0.0;

	for(k=N/2+1; k<N; k++){
		Diff[k]	= Dk*((double)(k - N));
	}

	if(contin == -1)
	{
		for(k=0; k<N; k++)
		{
			RkCur[k] = 0.0 + I*0.0;
			VkCur[k] = 0.0 + I*0.0;
		}
		RkCur[0] = 1.0 + I*0.0;
		VkCur[0] = 0.0 + I*0.0;
		
		read_etak();
		read_psik();
		find_z();
		find_R();
		find_V();
		
	}
	else{
		read_RV();
	}

/*	for (int j = 0; j < N; j++)
	{
		RuCur[j] = 0.1 * cexp(-I * k0 * Dk * j * du);
		VuCur[j] = -I * 0.1 * grav/(k0 * Dk) * cexp(-I * k0 * Dk * j * du);
	}

	fftw_execute(RuCur_2_RkCur);
	fftw_execute(VuCur_2_VkCur);

	for (int j = 0; j < N; j++)
	{
		RkCur[j] = RkCur[j]/N;
		VkCur[j] = VkCur[j]/N;
	}*/

	printf("creal(RkCur[0]) = %e\n", creal(RkCur[0]));

	for(k=0; k<N; k++){
		RkOld[k] = RkCur[k];
		VkOld[k] = VkCur[k];
		RkNew[k] = 0.0 + I*0.0;
		VkNew[k] = 0.0 + I*0.0;
	}
	
}
void open_files(void){
	if(contin==0){
		energy = fopen("./conserv/energy.dat","w");	
	}
	else{
		energy = fopen("./conserv/energy.dat","a");
	}
}


void recover_z(void){
	int i, j;
	fftw_complex sum;

	dZk[0] = 1.0;
	
	for(i=1; i<N/2; i++){
		sum = 0.0 + I*0.0;
		for(j=1; j<i; j++){
			sum = sum + RkCur[N-j]*dZk[N-(i-j)];
		}
		dZk[N-i] = - RkCur[N-i] - sum;
	}	
	for(i=1; i<=N/2; i++){
		dZk[i] = 0.0 + I*0.0;
	}
	fftw_execute(dZk_2_dZu_);

	
	for(i=1; i<N/2; i++){
		dXTk[N-i] = dZk[N-i]/2.0;
		XTk[N-i] = I*dXTk[N-i]/(((double)i)*Dk);
		Yk[N-i] = dZk[N-i]/2.0/(((double)i)*Dk);
		dXTk[i] = conj(dXTk[N-i]);
		XTk[i] = conj(XTk[N-i]);
		Yk[i] = conj(Yk[N-i]);
	}
	dXTk[0] = 0.0 + I*0.0;
	dXTk[N/2] = 0.0 + I*0.0;
	XTk[0] = 0.0 + I*0.0;
	XTk[N/2] = 0.0 + I*0.0;
	Yk[0] = 0.0 + I*0.0;
	Yk[N/2] = 0.0 + I*0.0;


	fftw_execute(dXTk_2_dXTu_);
	fftw_execute(XTk_2_XTu_);

	sum = 0.0 + I*0.0;
	for(i=1; i<N/2; i++){
		sum = sum + ((double)i)*Dk*Yk[i]*conj(Yk[i]);
	}
	Yk[0] = - 2.0*sum;
	

	fftw_execute(Yk_2_Yu_);

	EtaMax = 0.0;
	EtaMin = 0.0;


	for(i = 0; i < N; i++){
		if( EtaMax < creal(Yu[i]) ){
			EtaMax = creal(Yu[i]);
		}

		if( EtaMin > creal(Yu[i]) ){
			EtaMin = creal(Yu[i]);
		}
	}

}

void recover_psi(void){
	int u, i;
	
	fftw_execute(VkCur_2_VuCur_);

	for(u = 0; u < N; u++){
		dPhiu[u] = -I * VuCur[u] * dZu[u];
		dPsiu[u] = creal(dPhiu[u]);
	}
	fftw_execute(dPhiu_2_dPhik);
	fftw_execute(dPsiu_2_dPsik);

	for(i = 0; i < N; i++){
		dPhik[i] = dPhik[i] * norm;
		dPsik[i] = dPsik[i] * norm;
	}

// ========================= Distribution ==============
	HdPsik[0] = 0.0 + I * 0.0;
	for(i = 1; i < N; i++){
		HdPsik[i] = I * dPsik[i];
		HdPsik[N-i] = -I * dPsik[N-i];
	}
	HdPsik[N/2] = 0.0 + I * 0.0;
	fftw_execute(HdPsik_2_HdPsiu_);
// ========================= Distribution ==============
	for(i = 1; i < N/2; i++){
		Phik[i] = - I * dPhik[i] / ( ((double)i) * Dk );
		Psik[i] = - I * dPsik[i] / ( ((double)i) * Dk );
		Phik[N-i] = I * dPhik[N-i] / ( ((double)i) * Dk );
		Psik[N-i] = I * dPsik[N-i] / ( ((double)i) * Dk );
	}
	Phik[0] = 0.0 + I * 0.0;
	Phik[N/2] = 0.0 + I * 0.0;
	Psik[0] = 0.0 + I * 0.0;
	Psik[N/2] = 0.0 + I * 0.0;
	fftw_execute(Phik_2_Phiu_);
	fftw_execute(Psik_2_Psiu_);


}



// void recover_psi(void){
// 	int u, i;
// 	for(u=0; u<N; u++){
// 		dPhiu[u] = -I*VuCur[u]*dZu[u];
// 	}
// 	fftw_execute(dPhiu_2_dPhik);

// 	for(i=1; i<N/2; i++){
// 		dPhik[N-i] = dPhik[N-i]*norm;
// 		dPhik[i] = 0.0 + I*0.0;
// 		dPsik[N-i] = dPhik[N-i]/2.0;
// 		dPsik[i] = conj(dPsik[N-i]);
// 		Psik[N-i] = I*dPsik[N-i]/((double)i);
// 		Psik[i] = conj(Psik[N-i]);
// 	}
// 	dPhik[N/2] = 0.0 + I*0.0;
// 	dPhik[0] = 0.0 + I*0.0;
// 	dPsik[N/2] = 0.0 + I*0.0;
// 	dPsik[0] = 0.0 + I*0.0;
// 	Psik[N/2] = 0.0 + I*0.0;
// 	Psik[0] = 0.0 + I*0.0;	

// }

void calc_energy(void){
	int u;
	
	Epot = 0.0;
	Ekin = 0.0;
// ========================= Distribution ==============	
	EpotMax = 0.0;
	EpotMin = 0.0;
	EkinMax = 0.0;
	EkinMin = 0.0;
// ========================= Distribution ==============	
	

	for(u = 0; u < N; u++ ){
		Epot += grav/2.0 * pow( creal(Yu[u]), 2.0) * ( 1.0 + creal(dXTu[u]) );
		Ekin += 1.0/2.0 * creal( Psiu[u] ) * creal( HdPsiu[u] );
		
		Epotu[u] = grav/2.0 * pow( creal(Yu[u]), 2.0) * ( 1.0 + creal(dXTu[u]) ) / grav;
		Ekinu[u] = 1.0/2.0 * creal( Psiu[u] ) * creal( HdPsiu[u] ) / grav;


		if( EpotMax < Epotu[u]){
			EpotMax = Epotu[u];
		}
		if( EpotMin > Epotu[u]){
			EpotMin = Epotu[u];
		}
		if( EkinMax < Ekinu[u]){
			EkinMax = Ekinu[u];
		}
		if( EkinMin > Ekinu[u]){
			EkinMin = Ekinu[u];
		}

	}
	
	Epot = Epot / ( (double)N ) / grav;
	Ekin = Ekin / ( (double)N ) / grav;
	
	Etot = Epot + Ekin;

}
void calc_energy_old(void){
	int u;
	
	Epot = 0.0;
	Ekin = 0.0;
	
	fftw_complex sum = 0.0 + I * 0.0;

	for(u = 0; u < N; u++ ){
		Epot += pow( creal(Yu[u]), 2.0) * ( 1.0 + creal(dXTu[u]) );
	}
	for(u = 1; u < N/2; u++){
		sum += ( (double)u ) * Dk * Psik[u] * conj( Psik[u] );
	}
	sum = 2.0 * sum;
	
	Ekin = 1.0/2.0 * creal(sum) / grav;
	Epot = Epot * grav / 2.0 / ( (double)N ) / grav;
	// Ekin = 1.0/2.0*creal(sum)*Length/grav; //    /grav !!!
	// Epot = Epot*grav/2.0*Length/((double)N)/grav; //    /grav !!!
	Etot = Epot + Ekin;

}
void calc_momentum(void){
	int u;
	Momx = 0.0;
	Momy = 0.0;

// ========================= Distribution ==============
	MomxMax = 0.0;
	MomxMin = 0.0;
	MomyMax = 0.0;
	MomyMin = 0.0;
// ========================= Distribution ==============
	for(u = 0; u < N; u++){
		
		Momx += creal( Psiu[u] ) * cimag( dZu[u] );
		Momy += creal( Psiu[u] ) * creal( dZu[u] );
		
		Momxu[u] = creal( Psiu[u] ) * cimag( dZu[u] );
		Momyu[u] = creal( Psiu[u] ) * creal( dZu[u] );

		if( MomxMax < Momxu[u]){
			MomxMax = Momxu[u];
		}
		if( MomxMin > Momxu[u]){
			MomxMin = Momxu[u];
		}
		if( MomyMax < Momyu[u]){
			MomyMax = Momyu[u];
		}
		if( MomyMin > Momyu[u]){
			MomyMin = Momyu[u];
		}

		// Momx += creal(Psiu[u])*cimag(1.0/RuCur[u]);
		// Momy += creal(Phiu[u])*creal(1.0/RuCur[u]);
		// Momx += creal(Phiu[u])*cimag(1.0/RuCur[u]);
		// printf("%e\n", creal(Phiu[u])*creal(1.0/RuCur[u]) );
	}
	// Momx = Momx*2.0*M_PI/((double)N);
	// Momy = Momy*2.0*M_PI/((double)N);
	// Momx = Momx*Length/((double)N);
	// Momy = Momy*Length/((double)N);
	Momx = Momx / ( (double)N );
	Momy = Momy / ( (double)N );
}
void write_conserv(){
	fprintf(energy, "%e   %.18e   %.18e   %.18e   %.18e   %.18e\n", Time, Etot, Epot, Ekin, Momx, Momy );
}
void close_files(void){
	fclose(energy);
}


double calc_position(void){
	int i;
	double VMaxCur = 0.0;
	double xMaxCur = 0.0;

	fftw_execute(VkCur_2_VuCur_);
	for(i = 0; i < N; i++){
		if( sqrt( creal( VuCur[i] * conj(VuCur[i]) ) ) > VMaxCur ){
			VMaxCur = sqrt( creal( VuCur[i] * conj(VuCur[i]) ) );
			xMaxCur = du * i;
		}
	}

	// flag_gr_3 = 0;
	// if(VMaxCur >= Max_thr){
	// 	flag_gr_3 = 1;
	// 	ncalc_gr_3 = ncalc_gr_3 + 1;
	// }


	// printf("xMaxCur = %e\n", xMaxCur);
	return xMaxCur;


}


void correct_velocity(void){
	if(first_time == 1){
		first_time = 0.0;
		// UAvTp = (PosAvTpNew - x0) / Tp;
		// Vgr = Vgr + UAvTp;
	}
	else{

		UAvTp = (PosAvTpNew - PosAvTpOld) / Tp;
		Vgr = Vgr + UAvTp;

	}

	printf("PosAvTpNew = %e   PosAvTpOld = %e   UAvTp = %e\n", PosAvTpNew, PosAvTpOld, UAvTp);
	printf("Vgr = %e   Vgr0 = %e\n", Vgr, 0.5 * sqrt(grav / (Dk * kx0) ));

	PosAvTpOld = PosAvTpNew;
	PosAvTpNew = 0.0;
	
	
}

void damping_in_x_space(void){
	int i, k;
	double u;
	double x;

	fftw_execute(RkCur_2_RuCur_);
	fftw_execute(VkCur_2_VuCur_);

	// if( first_time == 0){
	// 	dx_correction = PosAvTpOld - Length / 2.0;
	// }
	

	for(i = 0; i < N; i++){
		x = du * i;// - dx_correction;
		u = x / Length * 2.0 * M_PI;
		RuCur[i] = 1.0 + (RuCur[i]-1.0) * exp( - Dcoeff * pow( cos(u/2.0), 6.0 ) * tau );

		VuCur[i] = VuCur[i] * exp( - Dcoeff * pow( cos(u/2.0), 6.0 ) * tau );

		// RuCur[i] = RuCur[i] * exp( - Dcoeff * pow( cos(u/2.0), 6.0 ) * tau );
		// VuCur[i] = VuCur[i] * exp( - Dcoeff * pow( cos(u/2.0), 6.0 ) * tau );

		
	}

	
	fftw_execute(RuCur_2_RkCur);
	fftw_execute(VuCur_2_VkCur);

	for(i = 0; i < N; i++){
		RkCur[i] = RkCur[i] / ( (double)N );
		VkCur[i] = VkCur[i] / ( (double)N );
	}

	
	RkCur[0] = 1.0;
	VkCur[0] = 0.0 + I * 0.0;
	
		
	for(k = N/2 + 1; k < N; k++){
		RkNew[k] = RkCur[k];
		RkOld[k] = RkCur[k];
		
		VkNew[k] = VkCur[k];
		VkOld[k] = VkCur[k];
	}
	for(k = 1; k <= N/2; k++){
		RkOld[k] = 0.0 + I * 0.0;
		RkCur[k] = 0.0 + I * 0.0;
		RkNew[k] = 0.0 + I * 0.0;

		VkOld[k] = 0.0 + I * 0.0;
		VkCur[k] = 0.0 + I * 0.0;
		VkNew[k] = 0.0 + I * 0.0;
	}

}


int main(void){
	int k, u;
	
	create_plan();
	
	initial();

	tCur = clock();

	//open_files();

	recover_z();
	recover_psi();

	//calc_energy();
	//calc_momentum();
	//write_conserv();
	//calc_steepness();
	write_evolution();
	write_restart();
	Right(rhsR1, rhsV1);

	while(Time <= FinTime)
	{
		
		Right(rhsR1, rhsV1);
		Shift(rhsR1, rhsV1, 0.5);
		Right(rhsR2, rhsV2);
		Shift(rhsR2, rhsV2, 0.5);
		Right(rhsR3, rhsV3);
		Shift(rhsR3, rhsV3, 1.0);
		Right(rhsR4, rhsV4);
		
		
		for(k = N/2 + 1; k < N; k++){
			RkNew[k] = RkOld[k] + tau * ( rhsR1[k] + 2.0 * rhsR2[k] + 2.0 * rhsR3[k] + rhsR4[k] ) / 6.0;
			VkNew[k] = VkOld[k] + tau * ( rhsV1[k] + 2.0 * rhsV2[k] + 2.0 * rhsV3[k] + rhsV4[k] ) / 6.0;
		}	

		RkNew[0] = 1.0;
		RkOld[0] = 1.0;
		RkCur[0] = 1.0;
		VkNew[0] = 0.0 + I * 0.0;
		VkOld[0] = 0.0 + I * 0.0;
		VkCur[0] = 0.0 + I * 0.0;
		
		for(k = N/2 + 1; k < N; k++){
			RkOld[k] = RkNew[k];
			RkCur[k] = RkNew[k];
			VkOld[k] = VkNew[k];
			VkCur[k] = VkNew[k];
		}
		for(k = 1; k <= N/2; k++){
			RkOld[k] = 0.0 + I * 0.0;
			RkCur[k] = 0.0 + I * 0.0;
			RkNew[k] = 0.0 + I * 0.0;
			VkOld[k] = 0.0 + I * 0.0;
			VkCur[k] = 0.0 + I * 0.0;
			VkNew[k] = 0.0 + I * 0.0;
		}

/*		if(Time < 200000.0)
		{
		 	damping_in_x_space();
		}*/
		


		Time += tau;
		nstep += 1;

/*		if( !( nstep % freq_ncalc_pos ) ){
			PosCur = calc_position();
			nstep_calc_pos = nstep_calc_pos + 1;
			if(nstep_calc_pos < Period_calc_pos){
				// if(flag_gr_3 == 1){
				// 	PosAvTpNew = PosAvTpNew + PosCur;	
				// }
				PosAvTpNew = PosAvTpNew + PosCur / Period_calc_pos;
				
			}
			else{
				// if(flag_gr_3 == 1){
					// PosAvTpNew = PosAvTpNew + PosCur;	
				// }
				// PosAvTpNew = PosAvTpNew / ncalc_gr_3;
				PosAvTpNew = PosAvTpNew + PosCur / Period_calc_pos;

				correct_velocity();
				
				
				ncalc_gr_3 = 0;
				nstep_calc_pos = 0;
			}
		}*/

/*		if( ! ( nstep % print_Enrg ) ){
			recover_z();
			recover_psi();
			//calc_energy();
			//calc_momentum();
			//write_conserv();
			//calc_steepness();

			nstep_Ev += 1;
			write_evolution();

		} */

		if(!(nstep%print_Restart))
		{
			recover_z();
			recover_psi();
			//calc_energy();
			//calc_momentum();
			//write_conserv();
			//calc_steepness();

			nstep_Ev += 1;
			nstep_Restart +=1;
			write_restart();
			write_evolution();
			printf("Time = %e   FinTime = %e   dEtaMax = %f\n", Time, FinTime, dEtaMax);
			
		}
		
		
	
		
		

		
	};
	
	
	recover_z();
	recover_psi();
	//calc_energy();
	//calc_momentum();
	//write_conserv();
	//calc_steepness();
	
	nstep_Ev += 1;
	write_evolution();
	
	nstep_Restart +=1;
	write_restart();
	printf("Time = %f   FinTime = %f   nstep_Restart = %d   nstep_Ev = %d\n", Time, FinTime, nstep_Restart, nstep_Ev);

	destroy_plan();
	// close_files();
	tCur = clock()-tCur;
	printf ("It took me %d clicks (%f seconds).\n",(int)tCur,((float)tCur)/CLOCKS_PER_SEC);
	return 0;
}

void write_evolution(void){
	int k;
	double uCur;
	char file_name1[40];
	char file_name2[40];
	// char file_name2[40];
	// char file_name3[40];
	char file_name4[40];
	char file_name5[40];
	char file_name6[40];

	char Name[50];
	
	sprintf(file_name1, "./R/%d_k.dat", nstep_Ev);
	pfileR = fopen(file_name1,"w");

	// sprintf(file_name2, "./en_mom/%d_x.dat", nstep_Ev);
	// pfileEnMom = fopen(file_name2,"w");

	sprintf(file_name5, "./surf/%d_x.dat", nstep_Ev);
	surf = fopen(file_name5,"w");

	sprintf(file_name6, "./psi/%d_x.dat", nstep_Ev);
	psi_data = fopen(file_name6,"w");

	// sprintf(file_name3, "./Ru/%d_x.dat", nstep_Ev);
	// data_Ru = fopen(file_name3,"w");

	//sprintf(file_name4, "./Vu/%d_x.dat", nstep_Ev);
	//data_Vu = fopen(file_name4,"w");
	

	fftw_execute(RkCur_2_RuCur_);
	fftw_execute(VkCur_2_VuCur_);
	
	fprintf (pfileR, "%e\t%e\n", 0.0, sqrt(pow(creal(RkCur[0]), 2.0) + pow(cimag(RkCur[0]), 2.0))  );

	for(k = 1; k < N/2; k++){
		fprintf (pfileR, "%e\t%e\n", -k*Dk, sqrt(pow(creal(RkCur[N-k]), 2.0) + pow(cimag(RkCur[N-k]), 2.0)) );
	}

	for(k = 0; k < N; k++){
		uCur = ( (double)k ) * du;
		//fprintf (surf, "%20.12e\t%20.12e\t%20.12e\n", uCur + creal(XTu[k]), creal(Yu[k]), creal(Steepness_u[k]) );
		fprintf (surf, "%20.12e\t%20.12e\t%20.12e\n", uCur + creal(XTu[k]), creal(Yu[k]), creal(Phiu[k]));

		// fprintf (data_Ru, "%20.12e   %20.12e   %20.12e\n", uCur + creal(XTu[k]), sqrt(creal( RuCur[k] * conj(RuCur[k]) ) ) , creal(RuCur[k]));
		//fprintf (data_Vu, "%20.12e   %20.12e   %20.12e\n", uCur + creal(XTu[k]), sqrt(creal( VuCur[k] * conj(VuCur[k]) ) ) , creal(VuCur[k]));
		fprintf (psi_data, "%20.12e   %20.12e\n", uCur + creal(XTu[k]), creal(Psiu[k]));
		// fprintf (pfileEnMom, "%20.12e   %20.12e   %20.12e   %20.12e   %20.12e\n", uCur + creal(XTu[k]), Epotu[k], Ekinu[k], Momxu[k], Momyu[k] );
	}
	fclose(pfileR);
	fclose(surf);

	// fclose(data_Ru);
	// fclose(data_Vu);


	fclose(psi_data);
	// fclose(pfileEnMom);
	
}
void Right(fftw_complex *rhsR, fftw_complex *rhsV){
	int k, u;

	fftw_execute(RkCur_2_RuCur_);
	fftw_execute(VkCur_2_VuCur_);

	for(u = 0; u < N; u++){
		Uu[u] = VuCur[u] * conj(RuCur[u]) + conj(VuCur[u]) * RuCur[u];
		dVVcu[u] = VuCur[u] * conj(VuCur[u]);
	}

	fftw_execute(Uu_2_Uk);
	fftw_execute(dVVcu_2_dVVck);

	
	for(k = 1; k <= N/2; k++){
		dVVck[k] = 0.0 + I * 0.0;
		dUk[k] = 0.0 + I * 0.0;
		Uk[k] = 0.0 + I * 0.0;

		dRk[k] = 0.0 + I * 0.0;//I*Diff[k]*RkCur[k];
		dVk[k] = 0.0 + I * 0.0;//I*Diff[k]*VkCur[k];
	}	
	Uk[0] = Uk[0] * norm / 2.0;
	dVVck[0] = 0.0 + I * 0.0;
	dUk[0] = 0.0 + I * 0.0;
	dRk[0] = 0.0 + I * 0.0;
	dVk[0] = 0.0 + I * 0.0;
	for(k = N/2 + 1; k < N; k++){
		Uk[k] = Uk[k] * norm; 
		dUk[k] = I * Diff[k] * Uk[k]; // нормирована в предыдущей строке 
		dVVck[k] = I * Diff[k] * dVVck[k] * norm;
		
		dRk[k] = I * Diff[k] * RkCur[k];
		dVk[k] = I * Diff[k] * VkCur[k];
	}
	fftw_execute(Uk_2_Uu_);
	fftw_execute(dVVck_2_dVVcu_);
	fftw_execute(dUk_2_dUu_);
	fftw_execute(dRk_2_dRu_);
	fftw_execute(dVk_2_dVu_);

	for(u = 0; u < N; u++){
		rhsRu[u] = Vgr * dRu[u] + I * ( Uu[u] * dRu[u] - dUu[u] * RuCur[u] );
		rhsVu[u] = Vgr * dVu[u] + I * ( Uu[u] * dVu[u] - RuCur[u] * dVVcu[u] ) + grav * ( RuCur[u] - 1.0 );
	}

	// for(u=0; u<N; u++){
	// 	rhsRu[u] = I*(Uu[u]*dRu[u] - dUu[u]*RuCur[u]);
	// 	rhsVu[u] = I*(Uu[u]*dVu[u] - RuCur[u]*dVVcu[u]) + grav*(RuCur[u] - 1.0);
	// }

	fftw_execute(rhsRu_2_rhsRk);
	fftw_execute(rhsVu_2_rhsVk);

	for(k = 0; k <= N/2; k++){
		// rhsRk[k] = 0.0 + I*0.0;
		// rhsVk[k] = 0.0 + I*0.0;
		rhsR[k] = 0.0 + I *0.0;
		rhsV[k] = 0.0 + I *0.0;
	}
	for(k=N/2+1; k<N; k++){
		// rhsRk[k] = rhsRk[k]*norm;
		// rhsVk[k] = rhsVk[k]*norm;
		rhsR[k] = rhsRk[k]*norm;
		rhsV[k] = rhsVk[k]*norm;
		// printf("%d   %e   %e\n", k, creal(rhsVk[k]), cimag(rhsVk[k]));
	}

	// for(k=k0-GICdisp; k<=k0+GICdisp; k++){
 // 		rhsR[N-k] = rhsR[N-k] + alpha_p*exp(-pow((k-k0)*Dk, 2.0)/2.0/pow(GICdisp*Dk, 2.0))*RkCur[N-k];
 // 		rhsV[N-k] = rhsV[N-k] + alpha_p*exp(-pow((k-k0)*Dk, 2.0)/2.0/pow(GICdisp*Dk, 2.0))*VkCur[N-k];
 // 	}
	// printf("\n");
	// printf("\n");


		
}


void Shift(fftw_complex* rhsR, fftw_complex* rhsV, double coeffRK){
	int k;
	
	for(k=N/2+1; k<N; k++){
		RkCur[k] = RkOld[k] + tau*coeffRK*rhsR[k];
		VkCur[k] = VkOld[k] + tau*coeffRK*rhsV[k];
	}
}





void create_plan(void){
// 	RkCur_2_RuCur_ = fftw_plan_dft_1d(N, &RkCur[0], &RuCur[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	VkCur_2_VuCur_ = fftw_plan_dft_1d(N, &VkCur[0], &VuCur[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	RuCur_2_RkCur = fftw_plan_dft_1d(N, &RuCur[0], &RkCur[0], FFTW_FORWARD, FFTW_EXHAUSTIVE );
// 	VuCur_2_VkCur = fftw_plan_dft_1d(N, &VuCur[0], &VkCur[0], FFTW_FORWARD, FFTW_EXHAUSTIVE );

	
// 	Uu_2_Uk = fftw_plan_dft_1d(N, &Uu[0], &Uk[0], FFTW_FORWARD, FFTW_EXHAUSTIVE );
// 	Uk_2_Uu_ = fftw_plan_dft_1d(N, &Uk[0], &Uu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	dVVck_2_dVVcu_ = fftw_plan_dft_1d(N, &dVVck[0], &dVVcu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	dVVcu_2_dVVck = fftw_plan_dft_1d(N, &dVVcu[0], &dVVck[0], FFTW_FORWARD, FFTW_EXHAUSTIVE );

// 	dUk_2_dUu_ = fftw_plan_dft_1d(N, &dUk[0], &dUu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	dRk_2_dRu_ = fftw_plan_dft_1d(N, &dRk[0], &dRu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	dVk_2_dVu_ = fftw_plan_dft_1d(N, &dVk[0], &dVu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	rhsRu_2_rhsRk = fftw_plan_dft_1d(N, &rhsRu[0], &rhsRk[0], FFTW_FORWARD, FFTW_EXHAUSTIVE );

// 	rhsVu_2_rhsVk = fftw_plan_dft_1d(N, &rhsVu[0], &rhsVk[0], FFTW_FORWARD, FFTW_EXHAUSTIVE );
// 	dZk_2_dZu_ = fftw_plan_dft_1d(N, &dZk[0], &dZu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	dXTk_2_dXTu_ = fftw_plan_dft_1d(N, &dXTk[0], &dXTu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	XTk_2_XTu_ = fftw_plan_dft_1d(N, &XTk[0], &XTu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );

// 	Yk_2_Yu_ = fftw_plan_dft_1d(N, &Yk[0], &Yu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	dPsiu_2_dPsik = fftw_plan_dft_1d(N, &dPsiu[0], &dPsik[0], FFTW_FORWARD, FFTW_EXHAUSTIVE );
// 	Psik_2_Psiu_ = fftw_plan_dft_1d(N, &Psik[0], &Psiu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	dPhiu_2_dPhik = fftw_plan_dft_1d(N, &dPhiu[0], &dPhik[0], FFTW_FORWARD, FFTW_EXHAUSTIVE );

// 	Phik_2_Phiu_ = fftw_plan_dft_1d(N, &Phik[0], &Phiu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );
// 	Yu_2_Yk = fftw_plan_dft_1d(N, &Yu[0], &Yk[0], FFTW_FORWARD, FFTW_EXHAUSTIVE );


// 	Psiu_2_Psik = fftw_plan_dft_1d(N, &Psiu[0], &Psik[0], FFTW_FORWARD, FFTW_EXHAUSTIVE );
// 	dPhik_2_dPhiu_ = fftw_plan_dft_1d(N, &dPhik[0], &dPhiu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );	

// 	dYk_2_dYu_ = fftw_plan_dft_1d(N, &dYk[0], &dYu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );	

// // ========================= Distribution ==============
// 	HdPsik_2_HdPsiu_ = fftw_plan_dft_1d(N, &HdPsik[0], &HdPsiu[0], FFTW_BACKWARD, FFTW_EXHAUSTIVE );	
// // ========================= Distribution ==============

// =============================================================================================

	Etak_2_Etax = fftw_plan_dft_1d(N, &Etak[0], &Etax[0], FFTW_BACKWARD, FFTW_ESTIMATE);
	Etax_2_Etak = fftw_plan_dft_1d(N, &Etax[0], &Etak[0], FFTW_FORWARD, FFTW_ESTIMATE);

	Psix_2_Psik = fftw_plan_dft_1d(N, &Psix[0], &Psik[0], FFTW_FORWARD, FFTW_ESTIMATE);

	RkCur_2_RuCur_ = fftw_plan_dft_1d(N, &RkCur[0], &RuCur[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	VkCur_2_VuCur_ = fftw_plan_dft_1d(N, &VkCur[0], &VuCur[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	RuCur_2_RkCur = fftw_plan_dft_1d(N, &RuCur[0], &RkCur[0], FFTW_FORWARD, FFTW_ESTIMATE );
	VuCur_2_VkCur = fftw_plan_dft_1d(N, &VuCur[0], &VkCur[0], FFTW_FORWARD, FFTW_ESTIMATE );

	
	Uu_2_Uk = fftw_plan_dft_1d(N, &Uu[0], &Uk[0], FFTW_FORWARD, FFTW_ESTIMATE );
	Uk_2_Uu_ = fftw_plan_dft_1d(N, &Uk[0], &Uu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	dVVck_2_dVVcu_ = fftw_plan_dft_1d(N, &dVVck[0], &dVVcu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	dVVcu_2_dVVck = fftw_plan_dft_1d(N, &dVVcu[0], &dVVck[0], FFTW_FORWARD, FFTW_ESTIMATE );

	dUk_2_dUu_ = fftw_plan_dft_1d(N, &dUk[0], &dUu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	dRk_2_dRu_ = fftw_plan_dft_1d(N, &dRk[0], &dRu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	dVk_2_dVu_ = fftw_plan_dft_1d(N, &dVk[0], &dVu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	rhsRu_2_rhsRk = fftw_plan_dft_1d(N, &rhsRu[0], &rhsRk[0], FFTW_FORWARD, FFTW_ESTIMATE );

	rhsVu_2_rhsVk = fftw_plan_dft_1d(N, &rhsVu[0], &rhsVk[0], FFTW_FORWARD, FFTW_ESTIMATE );
	dZk_2_dZu_ = fftw_plan_dft_1d(N, &dZk[0], &dZu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	dXTk_2_dXTu_ = fftw_plan_dft_1d(N, &dXTk[0], &dXTu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	XTk_2_XTu_ = fftw_plan_dft_1d(N, &XTk[0], &XTu[0], FFTW_BACKWARD, FFTW_ESTIMATE );

	Yk_2_Yu_ = fftw_plan_dft_1d(N, &Yk[0], &Yu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	dPsiu_2_dPsik = fftw_plan_dft_1d(N, &dPsiu[0], &dPsik[0], FFTW_FORWARD, FFTW_ESTIMATE );
	Psik_2_Psiu_ = fftw_plan_dft_1d(N, &Psik[0], &Psiu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	dPhiu_2_dPhik = fftw_plan_dft_1d(N, &dPhiu[0], &dPhik[0], FFTW_FORWARD, FFTW_ESTIMATE );

	Phik_2_Phiu_ = fftw_plan_dft_1d(N, &Phik[0], &Phiu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	Yu_2_Yk = fftw_plan_dft_1d(N, &Yu[0], &Yk[0], FFTW_FORWARD, FFTW_ESTIMATE );


	Psiu_2_Psik = fftw_plan_dft_1d(N, &Psiu[0], &Psik[0], FFTW_FORWARD, FFTW_ESTIMATE );
	dPhik_2_dPhiu_ = fftw_plan_dft_1d(N, &dPhik[0], &dPhiu[0], FFTW_BACKWARD, FFTW_ESTIMATE );	

	dYk_2_dYu_ = fftw_plan_dft_1d(N, &dYk[0], &dYu[0], FFTW_BACKWARD, FFTW_ESTIMATE );	

// ========================= Distribution ==============
	HdPsik_2_HdPsiu_ = fftw_plan_dft_1d(N, &HdPsik[0], &HdPsiu[0], FFTW_BACKWARD, FFTW_ESTIMATE );	
// ========================= Distribution ==============

}

void destroy_plan(void){

	fftw_destroy_plan(Etak_2_Etax);
	fftw_destroy_plan(Etax_2_Etak);

	fftw_destroy_plan(Psix_2_Psik);

	fftw_destroy_plan(RkCur_2_RuCur_);
	fftw_destroy_plan(VkCur_2_VuCur_);
	fftw_destroy_plan(RuCur_2_RkCur);
	fftw_destroy_plan(VuCur_2_VkCur);

	fftw_destroy_plan(Uu_2_Uk);
	fftw_destroy_plan(Uk_2_Uu_);
	fftw_destroy_plan(dVVck_2_dVVcu_);
	fftw_destroy_plan(dVVcu_2_dVVck);

	fftw_destroy_plan(dUk_2_dUu_);
	fftw_destroy_plan(dRk_2_dRu_);
	fftw_destroy_plan(dVk_2_dVu_);
	fftw_destroy_plan(rhsRu_2_rhsRk);

	fftw_destroy_plan(rhsVu_2_rhsVk);
	fftw_destroy_plan(dZk_2_dZu_);
	fftw_destroy_plan(dXTk_2_dXTu_);
	fftw_destroy_plan(XTk_2_XTu_);

	fftw_destroy_plan(Yk_2_Yu_);
	fftw_destroy_plan(dPsiu_2_dPsik);
	fftw_destroy_plan(Psik_2_Psiu_);
	fftw_destroy_plan(dPhiu_2_dPhik);

	fftw_destroy_plan(Phik_2_Phiu_);
	fftw_destroy_plan(Yu_2_Yk);

	fftw_destroy_plan(Psiu_2_Psik);
	fftw_destroy_plan(dPhik_2_dPhiu_);

	fftw_destroy_plan(dYk_2_dYu_);

// ========================= Distribution ==============
	fftw_destroy_plan(HdPsik_2_HdPsiu_);
// ========================= Distribution ==============
}