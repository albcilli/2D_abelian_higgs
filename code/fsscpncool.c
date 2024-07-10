#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#define M_PI 3.1415926535897932

#define D 2
#define N 10

#include<limits.h>
#include<time.h>
#include"dSFMT.h"
#include"dSFMT.c"


//this is the status of the random number generator
dsfmt_t rng_status;

//random number in (0,1)
double myrand(void){
	double ris;
	ris=dsfmt_genrand_open_open(&rng_status);
	return ris;
	}

//initialize random generator
void initrand(unsigned int s){
	unsigned int seed=s;
	if(s==0){
		seed=((unsigned int) time(NULL)+10) % UINT_MAX;
		if(seed==0){
			seed=1;
			}
		}
	dsfmt_init_gen_rand(&rng_status, seed);
	return;
	}

double correlation(int n, int k, double m, double a[]){
	int i;
	double c;
	c = 0;
	for(i=0; i<(n-k); i++){
		c += (a[i] - m)*(a[i+k] - m);
		}
	c = c/(n-k);
	return c;
	}

double mean(int n, int term, int e, double a[]){
	int i;
	double m;
	m = 0;
	for(i=term; i<n; i++) m += pow(a[i], e);
	m = m/(n-term);
	return m;
	}

double blocking(int n, int term, int r, double m, double a[], int e){
	int i, j, nb;
	double mb, dm;
	nb = ((n-term)/r)*r;
	dm = 0;
	for(i=0; i<(nb/r); i++){
		mb = 0;
		for(j=0; j<r; j++){
			mb += pow(a[r*i+j + term], e);
			}
		dm += pow((mb/r)-m,2);
		}
	dm = dm*r*r/nb/nb;
	dm = sqrt(dm);
	return dm;
	}

double corrlen(double g0, double gp, double temp1, double temp2, int L){
	return (1/(2*sin(M_PI/L)))*sqrt(g0/gp - 1);
	}

double binder(double u2, double u2sq, double temp1, double temp2, int L){
	return u2sq/(u2*u2);
	}
	
double b2coeff(double q2, double q4, double temp1, double temp2, int L){
	return (-q4 + 3*q2*q2)/(12*q2);
	}

double chixisq(double qsq, double g0, double gp, double temp1, int L){
	return (qsq/(L*L))*(1/(4*sin(M_PI/L)*sin(M_PI/L)))*(g0/gp - 1);
	}
	
double fsschit(double qsqL, double qsqinf, double temp1, double temp2, int L){
	return (qsqL/qsqinf);
	}

double fsschim(double G0L, double G0inf, double temp1, double temp2, int L){
	return G0L/G0inf;
	}

double fssb2(double q2L, double q4L, double q2inf, double q4inf, int L){
	return ((-q4L + 3*q2L*q2L)/(12*q2L))*((12*q2inf)/(-q4inf + 3*q2inf*q2inf));
	}
	
double bootstrap(int n, int term, int r, int M, double a[], double b[], double c[], double d[], int nhelp, int e1, int e2, int e3, int e4, double(*func)(double, double, double, double, int), int L){
	int i, j, k, l, nb;
	double *help1, *help2, *help3, *help4, o1, o2, o3, o4, *mb, m, dm;
	nb = ((n-term)/r)*r;
	mb = malloc(sizeof(double)*M);
	help1 = malloc(sizeof(double)*nb);
	if(nhelp == 2) help2 = malloc(sizeof(double)*nb);
	if(nhelp == 3) {
		help2 = malloc(sizeof(double)*nb);
		help3 = malloc(sizeof(double)*nb);
		}
	if(nhelp == 4){
		help2 = malloc(sizeof(double)*nb);
                help3 = malloc(sizeof(double)*nb);
                help4 = malloc(sizeof(double)*nb);
		}
	m = 0;
	dm = 0;
	for(i=0; i<M; i++){
		o1 = 0;
		o2 = 0;
		o3 = 0;
		o4 = 0;
		for(j=0; j<(nb/r); j++){//check nb/r or nb/r+1
			l = (int)(myrand()*nb);
			for(k=0; k<r; k++){
				if((l + k) == nb) l = 0;		
				help1[r*j + k] = a[l + k + term];
				o1 += pow(help1[r*j + k], e1);
				if(nhelp == 2) {
					help2[r*j + k] = b[l + k + term];
					o2 += pow(help2[r*j + k], e2);
					}
				else if(nhelp == 3){
					help2[r*j + k] = b[l + k + term];
                                        o2 += pow(help2[r*j + k], e2);
					help3[r*j + k] = c[l + k + term];
                                        o3 += pow(help3[r*j + k], e3);
					}
				else if(nhelp == 4){
                                        help2[r*j + k] = b[l + k + term];
                                        o2 += pow(help2[r*j + k], e2);
                                        help3[r*j + k] = c[l + k + term];
                                        o3 += pow(help3[r*j + k], e3);
                                        help4[r*j + k] = d[l + k + term];
                                        o4 += pow(help4[r*j + k], e4);
					}
				else {
					o2 += pow(help1[r*j + k], e2);
					o3 += pow(help1[r*j + k], e3);
					o4 += pow(help1[r*j + k], e4);
					}
				//if(r*j + k == nb) break;
				}
			}
		o1 = o1/nb;
		o2 = o2/nb;
		o3 = o3/nb;
		o4 = o4/nb;
		mb[i] = func(o1, o2, o3, o4, L);
		m += mb[i];
		dm += mb[i]*mb[i];
		}
	m = m/M;
	dm = dm/M;
	dm = sqrt(dm - m*m);
	free(help1);
	free(mb);
	if(nhelp == 2) free(help2);
	if(nhelp == 3) {
		free(help2);
		free(help3);
		}
	if(nhelp == 4){
		free(help2);
                free(help3);
		free(help4);
                }
	return dm;
	}

int main(void){
	FILE *data, *datainf, *result, *resfss;
	int i, j, k, l, r, M, temp1, seed, nmeas, nwords, nmeasinf, nwordsinf, dim, diminf, L[D], Linf[D], term, err, nboot, Nc = 7;
	double *g0, *gp, *mu2, **qu, Chixisq[Nc], dChixisq[Nc], G0, dG0, Gp, Qu[Nc], dQu[Nc], Q2[Nc], Q4[Nc], Xi, dXi, Chi[Nc], dChi[Nc], B2[Nc], dB2[Nc], U, dU, Mu2, Mu2sq, beta, lambda, temp2, **ener, Ene[Nc], dE, *g0inf, **quinf, G0inf, Q2inf[Nc], Q4inf[Nc], G0fss, dG0fss, Chifss[Nc], dChifss[Nc], B2fss, dB2fss, **rplaq, RPlaq[Nc], IPlaq, dRPlaq[Nc], dIPlaq;
	char ifname[100], ofname[100], stemp[100], stempinf[100], fssname[100];
	err = scanf("%s\n", stemp);
	if(err != 1){//file path
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	//printf("%s\n", stemp);
	err = scanf("%s\n", stempinf);
        if(err != 1){//file path
                fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
	//printf("%s\n", stempinf);
	if(scanf("%d\n", &seed) != 1){//seed for randgen
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(scanf("%d	%d	%d	%d\n", &nmeas, &nwords, &nmeasinf, &nwordsinf) != 4){//rows and columns
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(scanf("%d	%d	%d\n", &r, &term, &M) != 3){ //blocking block size, termsteps, bootstrap resamplings
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(scanf("%d	%d	%lf	%lf\n", &L[0], &Linf[0], &beta, &lambda) != 4){//L, beta and lambda parameters
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if(((nmeas*24) != nwords)||((nmeasinf*24) != nwordsinf)){
		fprintf(stderr, "Problems in data file (cutted line), (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
	dim = L[0];
	diminf = Linf[0];
	for(i=1; i<D; i++) {
		L[i] = L[0];
		Linf[i] = Linf[0];
		dim *= L[i];
		diminf *= Linf[i];
		}
	//snprintf(ifname, sizeof(ifname), "s", stemp, L[0], beta, lambda);
	data = fopen(stemp, "r");
	if(data == NULL){
		fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	datainf = fopen(stempinf, "r");
        if(datainf == NULL){
                fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }

	snprintf(ofname, sizeof(ofname), "res_L_%.3f_%.3f.txt", beta, lambda);
	//snprintf(ofname, sizeof(ofname), "res_%d_beta_%.3f.txt", L[0], lambda);
	result = fopen(ofname, "a");
	if(result == NULL){
		fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	snprintf(fssname, sizeof(fssname), "fss_L_%.3f_%.3f.txt", beta, lambda);
	resfss = fopen(fssname, "a");
        if(resfss == NULL){
                fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
	initrand(seed);
	g0 = malloc(sizeof(double)*nmeas);
	gp = malloc(sizeof(double)*nmeas);
	mu2 = malloc(sizeof(double)*nmeas);
	qu = malloc(sizeof(double)*Nc);
	rplaq = malloc(sizeof(double)*Nc);
	ener = malloc(sizeof(double)*Nc);
	for(i=0; i<Nc; i++){
		qu[i] = malloc(sizeof(double)*nmeas);
		rplaq[i] = malloc(sizeof(double)*nmeas);
		ener[i] = malloc(sizeof(double)*nmeas);
		}
	//iplaq = malloc(sizeof(double)*nmeas);
	//qsq = malloc(sizeof(double)*nmeas);
	//ener = malloc(sizeof(double)*nmeas);
	for(i=0; i<nmeas; i++){
		err = fscanf(data, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", &g0[i], &gp[i], &mu2[i], &qu[0][i], &rplaq[0][i], &ener[0][i], &qu[1][i], &rplaq[1][i], &ener[1][i], &qu[2][i], &rplaq[2][i], &ener[2][i], &qu[3][i], &rplaq[3][i], &ener[3][i], &qu[4][i], &rplaq[4][i], &ener[4][i], &qu[5][i], &rplaq[5][i], &ener[5][i], &qu[6][i], &rplaq[6][i], &ener[6][i]);
		if(err != 24){ 
                	fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
                	exit(EXIT_FAILURE);
                	}
		//qsq[i] = (qu[i]*qu[i]);
		}
	
	g0inf = malloc(sizeof(double)*nmeasinf);
        quinf = malloc(sizeof(double)*Nc);
	for(i=0; i<Nc; i++) quinf[i] = malloc(sizeof(double)*nmeasinf);
	for(i=0; i< nmeasinf; i++){
		err = fscanf(datainf, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", &g0inf[i], &temp2, &temp2, &quinf[0][i], &temp2, &temp2, &quinf[1][i], &temp2, &temp2, &quinf[2][i], &temp2, &temp2, &quinf[3][i], &temp2, &temp2, &quinf[4][i], &temp2, &temp2, &quinf[5][i], &temp2, &temp2, &quinf[6][i], &temp2, &temp2);
		if(err != 24){
                        fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
                        exit(EXIT_FAILURE);
                        }
		//printf("%f	%f\n", g0inf[i], quinf[i]);
		}
	G0 = mean(nmeas, term, 1, g0);
	Gp = mean(nmeas, term, 1, gp);
	Mu2 = mean(nmeas, term, 1, mu2);
        Mu2sq = mean(nmeas, term, 2, mu2);
	for(i=0; i<Nc; i++){
		Ene[i] = mean(nmeas, term, 1, ener[i]);
		RPlaq[i] = mean(nmeas, term, 1, rplaq[i]);
		dRPlaq[i] = blocking(nmeas, term, r, RPlaq[i], rplaq[i], 1);
        	Qu[i] = mean(nmeas, term, 1, qu[i]);
        	Q2[i] = mean(nmeas, term, 2, qu[i]);
        	Q4[i] = mean(nmeas, term, 4, qu[i]);
        	Chi[i] = Q2[i]/dim;
		B2[i] = b2coeff(Q2[i], Q4[i], temp2, temp2, L[0]);
		}
	//IPlaq = mean(nmeas, term, 1, iplaq);
        //dIPlaq = blocking(nmeas, term, r, IPlaq, iplaq, 1);
	//E = mean(nmeas, term, 1, ener);
	//dE = blocking(nmeas, term, r, E, ener);
	

	dG0 = blocking(nmeas, term, r, G0, g0, 1);
	for(i=0; i<Nc; i++){
		dQu[i] = blocking(nmeas, term, r, Qu[i], qu[i], 1);
		dChi[i] = blocking(nmeas, term, r, Q2[i], qu[i], 2)/dim;
		dB2[i] = bootstrap(nmeas, term, r, M, qu[i], qu[i], qu[i], qu[i], 1, 2, 4, 0, 0, b2coeff, L[0]);
		}


	Xi = corrlen(G0, Gp, temp2, temp2, L[0]);
	U = binder(Mu2, Mu2sq, temp2, temp2, L[0]);
	//Chixisq = chixisq(Q2, G0, Gp, L[0]);
	
	dXi = bootstrap(nmeas, term, r, M, g0, gp, qu[0], qu[0], 2, 1, 1, 0, 0, corrlen, L[0]);
	dU = bootstrap(nmeas, term, r, M, mu2, mu2, qu[0], qu[0], 1, 1, 2, 0, 0, binder, L[0]);
	//dChixisq = bootstrap(nmeas, term, r, M, qsq, g0, gp, 3, 1, 1, 1, chixisq, L[0]);
	
	G0inf = mean(nmeasinf, term, 1, g0inf);
	for(i=0; i<Nc; i++){
		Q2inf[i] = mean(nmeasinf, term, 2, quinf[i]);
		Q4inf[i] = mean(nmeasinf, term, 4, quinf[i]);
		}
	//printf("%d	%d	%d	%f	%f\n", Linf[0], diminf, nmeasinf, G0, G0inf);	

	//B2fss = fssb2(Q2, Q4, Q2inf, Q4inf, L[0]);
	
	//printf("%f	%f	%f\n", G0fss, Chifss, B2fss);	
	nboot = nmeasinf;
	if(nmeasinf > nmeas) nboot = nmeas;
	G0fss = fsschim(G0, G0inf, temp2, temp2, L[0]);
	dG0fss = bootstrap(nboot, term, r, M, g0, g0inf, qu[0], qu[0], 2, 1, 1, 0, 0, fsschim, L[0]);
	//dB2fss = bootstrap(nboot, term, r, M, qu, qu, quinf, quinf, 4, 2, 4, 2, 4, fssb2, L[0]);
	
	for(i=0; i<Nc; i++){
		Chifss[i] = ((double)diminf/dim)*fsschit(Q2[i], Q2inf[i], temp2, temp2, L[0]);
		dChifss[i] = ((double)diminf/dim)*bootstrap(nboot, term, r, M,  qu[i], quinf[i], qu[i], qu[i], 2, 2, 2, 0, 0, fsschit, L[0]);
        	Chixisq[i] = chixisq(Q2[i], G0, Gp, temp2, L[0]);
        	dChixisq[i] = bootstrap(nmeas, term, r, M, qu[i], g0, gp, qu[i], 3, 2, 1, 1, 0, chixisq, L[0]);
		}
	
        //fprintf(result, "%d %.3f %.3f %f %f   %f %f   %f %f   %f %f   %.8f %.8f       %f %f   %f %f   %f %f   %f %f\n", L[0], beta, lambda, Xi, dXi, U, dU, Qu, dQu, G0, dG0, Chi, dChi, B2, dB2, Chixisq, dChixisq, RPlaq, dRPlaq, IPlaq, dIPlaq);
        //fprintf(resfss, "%d %.3f %.3f   %f %f   %f %f   %f %f   %f %f   %f %f   %f %f   %f %f\n", L[0], beta, lambda, Xi, dXi, U, dU, Qu, dQu, G0, dG0, Chi, dChi, B2, dB2, Chixisq, dChixisq);

	fprintf(result, "%d	%.2f	%.2f	%f	%f	%f	%f	%f	%f", L[0], beta, lambda, Xi, dXi, U, dU, G0, dG0);
	for(i=0; i<Nc; i++){
		//printf("%f	", Ene[i]);
		fprintf(result, "	%f	%f	%.8f	%.8f	%f	%f	%f	%f	%f	%f", Qu[i], dQu[i], Chi[i], dChi[i], B2[i], dB2[i], Chixisq[i], dChixisq[i], RPlaq[i], dRPlaq[i]);
	}
	//printf("\n");
	fprintf(result, "\n");
	fprintf(resfss, "%d	%.2f	%.2f	%f	%f	%.8f	%.8f", L[0], beta, lambda, Xi, dXi, G0fss, dG0fss);
	for(i=0; i<Nc; i++){
		fprintf(resfss, "	%.8f	%.8f", Chifss[i], dChifss[i]);
		}
	fprintf(resfss, "\n");
	//fprintf(result, "%d %.2f %.2f %.10f %.10f   %.10f %.10f   %.10f %.10f   %.10f %.10f	%.10f %.10f\n", L[0], beta, lambda, G0, dG0, Xi, dXi, Chi, dChi, Qu, dQu, E, dE); 
	free(g0);
	free(gp);
	free(mu2);
	for(i=0; i<Nc; i++){
		free(qu[i]);
		free(rplaq[i]);
		free(quinf[i]);
		free(ener[i]);
		}
	free(ener);
	free(qu);
	free(rplaq);
	//free(iplaq);
	//free(qsq);
	//free(ener);
	free(quinf);
	free(g0inf);
	fclose(data);
	fclose(result);
	fclose(resfss);
	return 0;
	}
