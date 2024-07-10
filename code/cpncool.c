#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#define M_PI 3.1415926535897932

#define N 10
#define D 2
#define IMPR 1

//#define DEBUG

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
	}

//calculating the module of a N component vector
double module(double complex vect[N]){
	int i;
	double mod;
	mod = 0;
	for(i=0; i<N; i++) mod += cabs(vect[i])*cabs(vect[i]);//calculating the vector square module
	return sqrt(mod);//returning the vector module
	}

//normalizing the vector with N components in position k
void normalize(int k, double complex **vect){
	int i;
	double mod;
	mod = 0;	
	for(i=0; i<N; i++) mod += cabs(vect[k][i])*cabs(vect[k][i]);//calculating the vector square module
	mod = sqrt(mod);
	for(i=0; i<N; i++) vect[k][i] = vect[k][i]/mod;//normalizing
	}

//definining coordinates in D dimension with lexicographic ordering, finding next and previous position
void geometry(int dim, int size[D], int **next, int **prev){
	int i, j, k, frac, temp, x[D], xp[D], xm[D];
	for(i=0; i<dim; i++){//the lexicographic indexing is defined as i = x_0 + x_1*L_0 + ... + L_0*...*L_{D-2}*x_{D-2}, with x_{D-1} = i/(L_0*...*L_{D-2})
		frac = 1;
		temp = 0;
		for(j=1; j<D+1; j++){//calculating the cartesian coordinates for each site "i" (i is the lexicographic index)
			frac *= size[D-j];//dividing factor
			x[D-j] = (i-temp)/(dim/frac);//calculate the D-j component subtracting terms in the rhs to i and dividing by the correct factor
			temp +=  x[D-j]*(dim/frac);//subtracting term to i in the previous definition of the lexicographic indexing
			xp[D-j] = x[D-j] + 1;//next position in D-j direction
			xm[D-j] = x[D-j] - 1;//previous position in D-j direction
			if(xp[D-j] == size[D-j]) xp[D-j] = 0;//PBC
			if(xm[D-j] < 0) xm[D-j] = size[D-j] - 1;//PBC
			}
		for(j=1; j<D+1; j++){//going back from cartesian coordinates to lexicographic indexing for next and prevoius position
			frac = 1;
			next[i][D-j] = 0;
			prev[i][D-j] = 0;
			for(k=1; k<D+1; k++){//using two for cycles because the j index defines the next and prev direction, while the k index runs on the components of each point
				frac *= size[D-k];//same dividing factor
				if(k == j) {//the coordinate k in direction j is the one from xp/xm if k==j, else is the same of the point x
					next[i][D-j] += (dim/frac)*xp[D-k];
					prev[i][D-j] += (dim/frac)*xm[D-k];
					continue;
					}
				next[i][D-j] += (dim/frac)*x[D-k];//summing all terms to get the lexicographic indexing: only the D-j coordinate is xp/xm while others remain x
				prev[i][D-j] += (dim/frac)*x[D-k];
				}	
			}
		}
	}

//initializing both gauge and scalar fields
void inizialize(int iflag, int dim, double *ds, double *dg, char *config, double complex **scalfld, double complex **gaugefld){
	int i, j;
	double a, b;
	FILE *ifp;
	if(iflag == 1){
		for(i=0; i<dim; i++){//hot starting
			for(j=0; j<D; j++) gaugefld[i][j] = cexp((1 - 2*myrand())*I);
			for(j=0; j<N; j++) scalfld[i][j] = (1 - 2*myrand()) + (1 - 2*myrand())*I;
			normalize(i, scalfld);	
			}
		}
	else if(iflag == 0){//cold starting
		for(i=0; i<dim; i++){
			for(j=0; j<D; j++) gaugefld[i][j] = 1;
			for(j=0; j<N; j++) scalfld[i][j] = 1 + 1*I;			
			normalize(i, scalfld);	
			}
		}
	else{//loading a previous configuration from file
		ifp = fopen(config, "r");
		if(ifp == NULL){
			fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		if(fscanf(ifp, "%lf	%lf\n", ds, dg) != 2){//loading metropolis parameters of the previous simulation
			fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		for(i=0; i<dim; i++){
			for(j=0; j<N; j++){//loading scalar field
				if(fscanf(ifp, "%lf+i*%lf\n", &a, &b) != 2){
					fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				scalfld[i][j] = a+ b*I;
				}
			for(j=0; j<D; j++){//loading gauge fields
				if(fscanf(ifp, "%lf+i*%lf\n", &a, &b) != 2){
					fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				gaugefld[i][j] = a + b*I;
				}
			}
		fclose(ifp);
		}
	}

//multiplicating two field components for a generic SU(2) matrix, j and k are the considered indexes
void su2multip(int j, int k, double delta, double complex field[N]){
	double complex temp[2];
	double a, b, c, d, ph1, ph2, mod;
	#ifdef DEBUG
		if(j==k || j>= N || k>=N){//error if the indexes are out of range
			fprintf(stderr, "Out of range, (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	#endif
	temp[0] = field[j];//saving actual field values
	temp[1] = field[k];
	a = 1 + (2.0*myrand()-1.0)*delta;//extracting four coefficients for Id and Pauli matrices, the su(2) matrix M is in the form a*Id + d*I*X - c*I*Y + b*I*Z
	b = (2.0*myrand()-1.0)*delta;//delta defines the distance from the Id matrix
	c = (2.0*myrand()-1.0)*delta;
	d = (2.0*myrand()-1.0)*delta;
	mod = sqrt(a*a + b*b + c*c + d*d);
	a /= mod;//normalizing the coefficient to have a su(2) matrix, where a^2 + b^2 + c^2 + d^2 = 1
	b /= mod;
	c /= mod;
	d /= mod;
	ph1 = (2.0*myrand() - 1.0)*delta;//random phase to have an u(2) matrix
	ph2 = (2.0*myrand() - 1.0)*delta;
	if(myrand() < 1.0/2.0){
		field[j] = (a + I*b)*temp[0] - (c - I*d)*temp[1];//su(2) matrix in the form [alpha -conj(beta)] with alpha = a+I*b, beta = c+I*d
		field[k] = (c + I*d)*temp[0] + (a - I*b)*temp[1];//			    [beta  conj(alpha)]
		field[j] *= cexp(I*ph1);//random phase
		field[k] *= cexp(I*ph2);
		}
	else{//inverse matrix for detailed balance
		field[j] = (a - I*b)*temp[0] + (c - I*d)*temp[1];
		field[k] = -(c + I*d)*temp[0] + (a + I*b)*temp[1];
		field[j] *= cexp(-I*ph1);
		field[k] *= cexp(-I*ph2);
		}
	}

//calculating the "flav" component of the force on scalar field
double complex sforce(int it, int flav, double beta, double complex **scalfld, double complex **gaugefld, int **next, int **prev){
	int j, ip, ipp, im, imm; 
	double complex force;
	double c1, c2;
	c1 = 1;//standard cp(N) action
	c2 = 0;
	if(IMPR == 1){//using  Symanzik improved action
		c1 = 4.0/3.0;
		c2 = -1.0/12.0;
		}
	force = 0;		

/*								.ipp(1)	
								^
								|
								|
								.ip(1)
								^	
								|
								|
    				       .imm(0) <--- .im(0) <--- .it ---> .ip(0) ---> .ipp(0)
								|
								|
								v
								.im(1)
								|
								|
								v
								.imm(1)				*/


													
	for(j=0; j<D; j++){//force for the scalar field is the sum in all direction of terms in the form U(x)*z(x+mu) and similars
		ip = next[it][j];//defining the next, prev, next to next and prev to prev positions on the lattice
		im = prev[it][j];
		ipp = next[next[it][j]][j];
		imm = prev[prev[it][j]][j];
		force += -2*beta*N*(c1*((scalfld[ip][flav])*(gaugefld[it][j]) + conj(gaugefld[im][j])*(scalfld[im][flav]))//considering the action in the form -2*N*beta*(c1Re(std) + c2Re(Sym)
		+ c2*((scalfld[ipp][flav])*(gaugefld[ip][j])*(gaugefld[it][j]) + conj(gaugefld[im][j])*conj(gaugefld[imm][j])*(scalfld[imm][flav])));//2 terms in std action(next+prev) +2 impr(2next+2prev)
		}//the action here is considered in the form S = -2betaN*Re(conj(z)*forcez)
	return force;
	}

//calculating the force on the "dir" gauge field
double complex gforce(int it, int dir, double beta, double lambda, double complex  **scalfld, double complex **gaugefld, int **next, int **prev){
	int j, k, ip, ipp, im, ipo, imo, idg;
	double complex force, staple;
	double c1, c2;
	c1 = 1;//standard cp(N) action
	c2 = 0;
	if(IMPR == 1){//using  Symanzik improved action
		c1 = 4.0/3.0;
		c2 = -1.0/12.0;
		}
	force = 0;
	ip = next[it][dir];
	im = prev[it][dir];
	ipp = next[next[it][dir]][dir];
	//imm = prev[prev[it][jt]][jt];
	/*			                                .ipp(1) 
                                                                ^
                                                                |
                                                                |
                                                                .ip(1)
                                                                ^       
                                                                | dir
                                                                |
                                       .imm(0) <--- .im(0) <--- .it ---> .ip(0) ---> .ipp(0)
                                                                |
                                                                |
                                                                v
                                                                .im(1)
                                                                |
                                                                |
                                                                v
                                                                .imm(1)                         */

	for(j=0; j<N; j++){//force on the gauge field is given by the sum on all flavours of terms in the form conj(fld(x+mu))*fld(x) and similars
		force += -2*beta*N*(c1*(scalfld[it][j])*conj(scalfld[ip][j]) + c2*((scalfld[it][j])*conj(scalfld[ipp][j])*conj(gaugefld[ip][dir]) //action in the form -2*N*beta*(c1Re(std) + c2Re(Sym)
		+ conj(scalfld[ip][j])*(scalfld[im][j])*conj(gaugefld[im][dir])));//one terms in std action(curr) and two in Symanzik(curr+prev)
		}//the action here is considered in the form S = -2betaN*Re(conj(U)*forceU)
	staple = 0;

		/*			  .ipo --->.
					  ^        |
					j |	   |
					  |  dir   v
					  .it <--- .ip
					  | 	   ^
					j |	   |
					  v        |    
					  .imo --->.idg		*/


	if(fabs(lambda) > 1/pow(10,10)){//calculating the plaquette term only if the lambda parameter is != 0
		for(k=1; k<D; k++){
			j = (dir+k)%D;//j runs on all directions ortogonal to dir
			ipo = next[it][j];//next and prev in ortogonal directions
			imo = prev[it][j];
			idg = next[imo][dir];//diagonal opposite point to it 
			staple += conj(gaugefld[ip][j])*gaugefld[ipo][dir]*gaugefld[it][j] + conj(gaugefld[imo][j])*gaugefld[imo][dir]*gaugefld[idg][j];//summing on all staples
			}
		force += -2*beta*lambda*staple;//action in the form S = conj(gauge*force), with force on the "dir" link = -2*beta*lambda*staple
		}
	return force;
	}

//calculating the system energy(action) on the given configuration
double energy(int dim, double beta, double lambda, double complex **scalfld, double complex **gaugefld, int **next){
	int i, j, k, ip, ipp, ipo;
	double c1, c2, ene, enetest;
	double complex staple;
	ene = 0;
	c1 = 1;//standard action
	c2 = 0;
	if(IMPR == 1){//improved action
		c1 = 4.0/3.0;
		c2 = -1.0/12.0;
		}
	for(i=0; i<dim; i++){
		for(j=0; j<D; j++){
			ip = next[i][j];//defining the next, prev, next to next and prev to prev positions on the lattice
			ipp = next[next[i][j]][j];
			for(k=0; k<N; k++){//calculating the energy(action) using the eq (11) in ref2
				ene += -2*N*beta*(c1*(creal(conj(gaugefld[i][j])*conj(scalfld[ip][k])*scalfld[i][k])) 
				+ c2*(creal(conj(gaugefld[ip][j])*conj(gaugefld[i][j])*conj(scalfld[ipp][k])*scalfld[i][k])));
				}
			ene += 2*N*beta*(c1 + c2);
			staple = 0;
			if(fabs(lambda) > 1/pow(10,10)) {//plaquette term
				for(k=j+1; k<D; k++){
					ipo = next[i][k];
					staple += conj(gaugefld[ip][k])*gaugefld[ipo][j]*gaugefld[i][k];
					}
				ene += -2*lambda*beta*(creal(conj(gaugefld[i][j])*staple));
				}
			}//overall action(energy) is in the form of eq (1) ref3
		}
	return ene;//returning the system energy
	}

//metroplis step for scalar field in the "it" site
int metros(int it, int dim, double ds, double beta, double lambda, double complex **scalfld, double complex **gaugefld, int **next, int **prev){
	int j, acc, temp0, temp1;
	double complex force[2], stemp[N], sold[N];
	double p, enetest;
	#ifdef DEBUG
		enetest = energy(dim, beta, lambda, scalfld, gaugefld, next);
	#endif
	acc = 0;//if acc = 1 the step is rejected, 0 if accepted
	for(j=0; j<N; j++) {
		stemp[j] = scalfld[it][j];//saving the actual scalar field configuration in position "it"
		sold[j] = scalfld[it][j];
		}
	temp0 = (int)(myrand()*N);//choosing two different indexes to multiplicate
	temp1 = (int)(myrand()*N);
	while(temp0 == temp1) temp1 = (int)(myrand()*N);
	su2multip(temp0, temp1, ds, stemp);//multipling two scalar field components with the su(2) defined above, ds defines how far it goes from identity
	for(j=0; j<N; j++) scalfld[it][j] = stemp[j];//saving the new configuration in the "it" site
	//metropolis test 
	#ifdef DEBUG
		enetest -= energy(dim, beta, lambda, scalfld, gaugefld, next);
	#endif
	p=0;
	force[0] = sforce(it, temp0, beta, scalfld, gaugefld, next, prev);//force here is calculated only in the two changed indexes
	force[1] = sforce(it, temp1, beta, scalfld, gaugefld, next, prev);//p = -S(test) -(-S(act)), with S = -2*beta*N(conj(field)*force), action from eq(11) ref2
	p = creal(-force[0]*conj(scalfld[it][temp0]) - force[1]*conj(scalfld[it][temp1]) + force[0]*conj(sold[temp0]) + force[1]*conj(sold[temp1]));//acceptance factor exponent
	#ifdef DEBUG
		if(fabs(p-enetest) > 1/pow(10, 10)) {//the acceptance exponent calculated in the two ways must be the same, if not the metro test fails
			fprintf(stderr, "Metro S energy test failed, (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	#endif
	if(myrand() < exp(p)) acc = 1;//acceptance test with exp(-(DeltaS)), if passed the configuration is updated
	else {
		for(j=0; j<N; j++) scalfld[it][j] = sold[j];//rejected move, keeping the old configuration
		}
	return acc;//returning 0 if the test fails and 1 if is successful
	}


//metroplis step for gauge field "dir" in the "it" site 
int metrog(int it, int dir, int dim, double dg, double beta, double lambda, double complex **scalfld, double complex **gaugefld, int **next, int **prev){
	int acc;
	double complex force, gtemp, gold;
	double p, enetest;
	acc = 0;//if acc = 1 the step is rejected, 0 if accepted
	#ifdef DEBUG
		enetest = energy(dim, beta, lambda, scalfld, gaugefld, next);
	#endif
	gold = gaugefld[it][dir];//saving the old configuration
	gtemp = cexp(dg*(1-2*myrand())*I)*gaugefld[it][dir];//adding a random phase in [-dg, dg]
	gaugefld[it][dir] = gtemp;//saving the new configuration
	//metropolis test
	#ifdef DEBUG
		enetest -= energy(dim, beta, lambda, scalfld, gaugefld, next); 
	#endif
	p=0; 
	force = gforce(it, dir, beta, lambda, scalfld, gaugefld, next, prev);//force here is calculated only in the dir link
	p = (creal(-force*conj(gaugefld[it][dir]) + force*conj(gold)));//p = -S(test) -(-S(act)), with S = -2*beta*N(conj(field)*force) -2*beta*lambda(conj(field)*staple), action from eq(1) ref3
	#ifdef DEBUG
		if(fabs(p-enetest) > 1/pow(10, 10)) {//the acceptance calculated in the two ways must be the same, if not the metro test fails
			fprintf(stderr, "Metro G energy test failed, (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	#endif
	if(myrand() < exp(p)) acc = 1; //acceptance test with exp(-(DeltaS)), if passed the configuration is updated
	else gaugefld[it][dir] = gold;//rejected move, keeping the old configuration
	return acc;//returning 0 if the test fails and 1 if is successful
	}

//microcanonical step for scalar field in position "it"
void micros(int it, int dim, double beta, double lambda, double complex **scalfld, double complex **gaugefld, int **next, int **prev){
	int i;
	double complex force[N], prod;
	double mod, microtest;
	#ifdef DEBUG
		microtest = energy(dim, beta, lambda, scalfld, gaugefld, next);//calculating the energy before the microcanonical step
	#endif
	for(i=0; i<N; i++) force[i] = sforce(it, i, beta, scalfld, gaugefld, next, prev);//force on the scalar field, writing the action in the form S = -2beta*N*Re(force*conj(fld))
	mod = 0;
	prod = 0;
	for(i=0; i<N; i++) {//calculating factors for eq (43a) ref1
		prod += conj(scalfld[it][i])*force[i];
		mod += conj(force[i])*force[i];
		}
	for(i=0; i<N; i++) scalfld[it][i] = 2*force[i]*creal(prod)/mod - scalfld[it][i];//reflection of scalar field according to eq (43a) and (44a) ref1
	#ifdef DEBUG
		microtest -= energy(dim, beta, lambda, scalfld, gaugefld, next);//calculating the energy difference after the step
		if(fabs(microtest) > 1/pow(10, 10)) {//microcanonical step preserves the energy, if the Delta E is not 0 the test fails
			fprintf(stderr, "Micro S energy test failed, (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	#endif
	}

//microcanonical step for gauge field "dir" in position "it"
void microg(int it, int dir, int dim, double beta, double lambda, double complex **scalfld, double complex **gaugefld, int **next, int **prev){
	double complex force;
	double microtest;
	//microcanonical test
	#ifdef DEBUG
		microtest = energy(dim, beta, lambda, scalfld, gaugefld, next);//calculating the energy before the microcanonical step
	#endif
	force = gforce(it, dir, beta, lambda, scalfld, gaugefld, next, prev);//force on the "dir" link, action in the form S = -2beta*N*Re(force*conj(fld)) -2*beta*lambda(conj(gauge)*staple)), eq(1) ref3
	gaugefld[it][dir] = 2*force*creal(conj(gaugefld[it][dir])*force)/(conj(force)*force) - gaugefld[it][dir];//reflection of field gauge as in eq (43a) and (44a) ref1
	#ifdef DEBUG
		microtest -= energy(dim, beta, lambda, scalfld, gaugefld, next);//calculating the energy difference after the step
		if(fabs(microtest) > 1/pow(10, 10)) {//microcanonical step preserves the energy, if the Delta E is not 0 the test fails
			fprintf(stderr, "Micro G energy test failed, in (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
	#endif
	}
/*
double energytest(int dim, double complex **scalfld, double complex **gaugefld, int **next){
        int i, j, k, ip, ipp;
        double c1, c2, ene;
        ene = 0;
        c1 = 1;//standard action
        c2 = 0;
        if(IMPR == 1){//improved action
                c1 = 4.0/3.0;
                c2 = -1.0/12.0;
                }
        for(i=0; i<dim; i++){
                for(j=0; j<D; j++){
                        ip = next[i][j];//defining the next, prev, next to next and prev to prev positions on the lattice
                        ipp = next[next[i][j]][j];
                        for(k=0; k<N; k++){//calculating the energy(action) using the eq (11) in ref2
                                ene += -2*(c1*(creal(conj(gaugefld[i][j])*conj(scalfld[ip][k])*scalfld[i][k]))
                                + c2*(creal(conj(gaugefld[ip][j])*conj(gaugefld[i][j])*conj(scalfld[ipp][k])*scalfld[i][k])));
                                }
			ene += 2*(c1 + c2);
                	}
		}
	ene = ene/(2*dim);
        return ene;//returning the system energy
        }
*/

//measuring considered observables
void measure(int dim, double beta, double lambda, int L[D], double complex **scalfld, double complex **gaugefld, int **next, FILE *ofp){
	double complex oper[N][N], poper[N][N], tempq, Plaq;
	double p[D], ft, g0, gp, mu2, Qu, qx, ener;
	int i, j, k, l, ip0, ip1, ip01, div, temp, x[D], err;	
	p[0] = 2*M_PI/L[0];//minimum momentum used in G(p) operator, eq(5) ref 3
	for(i=1; i<D; i++) p[i] = 0;//other components are 0
	for(j=0; j<N; j++){//initializing the operators sum(x)Qx and sum(x)Qpx
			for(k=0; k<N; k++){
					oper[j][k] = 0;
					poper[j][k] = 0;
				}
			}
	Qu = 0;
	Plaq = 0;	
	for(i=0; i<dim; i++){//calulating the sum on all sites of operators Qx = conj(field)*field - Id/N, and Qpx = exp(ipx)Qx
		ft = 0;
		div = 1;
		temp = 0;
		for(j=1; j<D+1; j++){//calulating cartesian coordinates of the point "i", same as geometry() function
			div *= L[D-j];
			x[D-j] = (i-temp)/(dim/div);
			temp += x[D-j]*(dim/div);
			//ft += x[D-l]*p[D-l];
			}
		ft += x[0]*p[0];//phase factor px in Qpx operator
		for(j=0; j<N; j++){
			for(k=0; k<N; k++){//calulating the operator sum(x)Qx and sum(x)Qpx
				oper[j][k] += conj(scalfld[i][j])*scalfld[i][k];
				poper[j][k] += cexp(ft*I)*conj(scalfld[i][j])*scalfld[i][k];
				if(j == k){//implementing the Id/N factor
					oper[j][k] -= 1.0/N;
					poper[j][k] -= cexp(ft*I)*(1.0/N);
					}
				}			
			}
		 if(D == 2){
                        ip0 = next[i][0];//next position in direction 0 
                        ip1 = next[i][1];//next position in direction 1
			tempq = clog((gaugefld[i][0])*(gaugefld[ip0][1])*conj(gaugefld[ip1][0])*conj(gaugefld[i][1]));
                        qx = cimag(tempq);//topological charge density calculation, eq (14) ref3
                        Qu += qx;//summing on all lattice to get the topological charge
			Plaq += (gaugefld[i][0])*(gaugefld[ip0][1])*conj(gaugefld[ip1][0])*conj(gaugefld[i][1]);
			}
		}
	Qu = Qu/(2*M_PI);//dividing factor
	Plaq = Plaq/dim;
	g0 = 0;//G0 (or susceptibility chi) is obtained as Tr(sum(x)Qx*sum(y)Qy) = sum(x)sum(y)Tr(QxQy) = sum(x)sum(y)G(x-y) = V*sum(x)G(x) = V*G0, see ref3 and handwritten notes
	gp = 0;//for Gpsame as above, with Qx -> exp(ipx)Qx
	mu2 = 0;//same as above but u2 is defined as sum(x)sum(y)Tr(QxQy), eq(6) ref3
	for(j=0; j<N; j++){
		for(k=0; k<N; k++){
			g0 += (oper[j][k]*oper[k][j]);
			gp += poper[j][k]*conj(poper[j][k]);//conj is needed on the 2nd operator in Gp
			}
		}
	mu2 = g0;//u2 is g0*dim
	g0 = g0/dim;//g0 is u2/dim
	gp = gp/dim;
	
	ener = energy(dim, beta, lambda, scalfld, gaugefld, next);//added for test
	
	err = fprintf(ofp, "%f	%f	%f	%f	%f	%f	", g0, gp, mu2, Qu, creal(Plaq), ener);//printing results on the file pointed by ofp
	if(err < 0){
		fprintf(stderr, "Error in file writing, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	}

void measurecool(int dim, double beta, double lambda, int L[D], double complex **scalfld, double complex **gaugefld, int **next, FILE *ofp, int end){
        double complex tempq, Plaq;
        double Qu, qx, ener;
        int i, ip0, ip1, ip01, err;
               
        Qu = 0;
        Plaq = 0;
	ener = energy(dim, beta, lambda, scalfld, gaugefld, next);//added for test
        if(D == 2){
		for(i=0; i<dim; i++){
        		ip0 = next[i][0];//next position in direction 0 
     	        	ip1 = next[i][1];//next position in direction 1
                	tempq = clog((gaugefld[i][0])*(gaugefld[ip0][1])*conj(gaugefld[ip1][0])*conj(gaugefld[i][1]));
                	qx = cimag(tempq);//topological charge density calculation, eq (14) ref3
                	Qu += qx;//summing on all lattice to get the topological charge
                	Plaq += (gaugefld[i][0])*(gaugefld[ip0][1])*conj(gaugefld[ip1][0])*conj(gaugefld[i][1]);
                	}
        	Qu = Qu/(2*M_PI);//dividing factor
        	Plaq = Plaq/dim;
		}
	if(end == 1){
        	err = fprintf(ofp, "%f	%f	%f\n", Qu, creal(Plaq), ener);//printing results on the file pointed by ofp
        	if(err < 0){
        		fprintf(stderr, "Error in file writing, (%s, %d)\n", __FILE__, __LINE__);
                	exit(EXIT_FAILURE);
			}
		}
	else{
		err = fprintf(ofp, "%f	%f	%f	", Qu, creal(Plaq), ener);//printing results on the file pointed by ofp
                if(err < 0){
                        fprintf(stderr, "Error in file writing, (%s, %d)\n", __FILE__, __LINE__);
                        exit(EXIT_FAILURE);
                        }
		}	
	}



//update function, the elementary step consist in a metropolis step on the full lattice, nmicr microcanonical steps on the full lattice and in the renormalization of the scalar field
void update(int dim, double ds, double dg, double beta, double lambda, double complex **scalfld, double complex **gaugefld, int **next, int **prev, double acc[2]){
	int i, j, dir, nmicr, accs, accg;
	nmicr = 5;// microcanonical vs metropolis steps ratio
	accs = 0;
	accg = 0;
	for(i=0; i<dim; i++){
		accs += metros(i, dim, ds, beta, lambda, scalfld, gaugefld, next, prev);//metropolis on both fields, with acceptance
		for(dir=0; dir<D; dir++) accg += metrog(i, dir, dim, dg, beta, lambda, scalfld, gaugefld, next, prev);
		}
	for(j=0; j<nmicr; j++){
		for(i=0; i<dim; i++){
			micros(i, dim, beta, lambda, scalfld, gaugefld, next, prev);//microcanonical steps on both fields
			for(dir=0; dir<D; dir++) microg(i, dir, dim, beta, lambda, scalfld, gaugefld, next, prev);
			}
		}
	for(i=0; i<dim; i++) normalize(i, scalfld);//normalizing the scalar field after the idec metropolis steps to ensure that is still normalized
	acc[0] = (double)accs/dim;
	acc[1] = (double)accg/(dim*D);
	}

void coolings(int it, double beta, double lambda, double complex **scalfld, double complex **gaugefld, int **next, int **prev){
	int j;
	double mod;
	double complex force[N];
        mod = 0;
        for(j=0; j<N; j++) {
		force[j] = sforce(it, j, beta, scalfld, gaugefld, next, prev);
                mod += conj(force[j])*force[j];
                }
	mod = sqrt(mod);
        for(j=0; j<N; j++) scalfld[it][j] = -force[j]/mod;

	}

void coolingg(int it, int dir, double beta, double lambda, double complex **scalfld, double complex **gaugefld, int **next, int **prev){
        double complex force; 

	force = gforce(it, dir, beta, lambda, scalfld, gaugefld, next, prev);
        gaugefld[it][dir] = -force/(sqrt(force*conj(force)));
        }

void updatecool(int dim, double beta, double lambda, double complex **scalfld, double complex **gaugefld, int **next, int **prev){
        int i, j, k, ncool;
        ncool = 5;
        for(k=0; k<ncool; k++){
                for(i=0; i<dim; i++){
                        coolings(i, beta, lambda, scalfld, gaugefld, next, prev);
                        for(j=0; j<D; j++) coolingg(i, j, beta, lambda, scalfld, gaugefld, next, prev);
                        }
                }
        for(i=0; i<dim; i++) normalize(i, scalfld);
        }



//saving the actual configuration and the metropolis parameters in the ofp file
void saveconfig(int dim, double ds, double dg, double complex **scalfld, double complex **gaugefld, FILE *ofp){
	int j, k, err;
	err = fprintf(ofp, "%f	%f\n", ds, dg);//saving metropolis paramteres
	if(err < 0){
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	for(k=0; k<dim; k++){	
		for(j=0; j<N; j++) {//saving scalar fields
			err = fprintf(ofp, "%.16f+i*%.16f\n", creal(scalfld[k][j]), cimag(scalfld[k][j]));//configurations saving with 16 decimals
			if(err < 0){
				fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			}
		for(j=0; j<D; j++) {//saving gauge fields
			err = fprintf(ofp, "%.16f+i*%.16f\n", creal(gaugefld[k][j]), cimag(gaugefld[k][j]));
			if(err < 0){
				fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			}
		}
	}

void copyconf(int dim, double complex **oldscal, double complex **oldgauge, double complex **newscal, double complex **newgauge){
	int j, k;
	for(k=0; k<dim; k++){
		for(j=0; j<N; j++){
			newscal[k][j] = oldscal[k][j];
			}
		for(j=0; j<D; j++){
			newgauge[k][j] = oldgauge[k][j];
			}
		}
	}

int main(void){
	int i, j, k, err, parity, temp, confn, seed, dim, nmeas, iflag, iterm, idec, L[D], **np, **nm;//defining all needed parameters and structures
	double complex **field, **gauge, microtests[N], microtestg[D], **coolfield, **coolgauge;//pointer for dim*N and dim*D arrays for scalar field and gauge field
	double beta, lambda, ds, dg, acc[2], accs, accg;//parameters for metropolis of both fields and action parameters
	FILE *configA, *configB, *configF, *log, *dati;//pointers to saving configurations file(config), saving odd, even and final configurations of 10, data and logs file
	char sdati[50], sconfigA[50], sconfigB[50], sconfigF[50], slog[50], stemp[50];
	/*if(fscanf(param, "%d	%d	%d	%d	%d	%lf	%lf	%lf	%lf\n", &seed, &nmeas, &iflag, &iterm, &idec, &beta, &lambda, &ds, &dg) != 9){
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}	*/
	err = scanf("%d\n", &seed);
	if(err != 1){//seed for randgen
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	err = scanf("%d	%d	%d	%d\n", &nmeas, &iterm, &idec, &iflag);
	if(err != 4){ //nmeas, termsteps, decorrsteps, cold(0) hot(1) or file(else) starting
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	err = scanf("%lf	%lf\n", &beta, &lambda);
	if(err != 2){//beta and lambda action parameters
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	err = scanf("%lf	%lf\n", &ds, &dg);
	if(err != 2){//scalar and gauge metropolis parameters
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}

	dim = 1;
	for(i=0; i<D; i++) {//lattice size in each direction, usually sqared
		err = scanf("%d\n", &L[i]);
		if(err != 1){
			fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}	
		dim *= L[i];//total lattice size
		}
	if(L[0] != L[1]){//squared lattice test
		fprintf(stderr, "Lattice not squared, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	err = scanf("%s\n", stemp);//string with path for config files
	if(err != 1){
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	snprintf(sconfigA, sizeof(sconfigA), "%s_%d_%.3f_%.3fA.dat", stemp, L[0], beta, lambda);
	snprintf(sconfigB, sizeof(sconfigB), "%s_%d_%.3f_%.3fB.dat", stemp, L[0], beta, lambda);
	snprintf(sconfigF, sizeof(sconfigF), "%s_%d_%.3f_%.3fF.dat", stemp, L[0], beta, lambda);
	err = scanf("%s\n", stemp);//string with path for data file
	if(err != 1){
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	snprintf(sdati, sizeof(sdati), "%s_%d_%.3f_%.3f.dat", stemp, L[0], beta, lambda);
	err = scanf("%s\n", stemp);//string with path for log file
	if(err != 1){
		fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	snprintf(slog, sizeof(slog), "%s_%d_%.3f_%.3f.dat", stemp, L[0], beta, lambda);
	confn = 99;
	if((iflag != 0) && (iflag != 1)){//using loaded configuration
		log = fopen(slog, "r");//reading in log file which was the last saved when using loaded configuration
		if(log == NULL){
			fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		err = fscanf(log, "Configuration = %d of 10\n", &confn); //reading in log file the number of the saved configuration of the 10 
		if(err != 1){
			fprintf(stderr, "Error in file reading, (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		fclose(log);
		}
	field = malloc(sizeof(double complex *)*dim);//structures allocation 
	gauge = malloc(sizeof(double complex *)*dim);
	np = malloc(sizeof(int *)*dim);
	nm = malloc(sizeof(int *)*dim);
	for(i=0; i<dim; i++){
		field[i] = malloc(sizeof(double complex)*N);
		gauge[i] = malloc(sizeof(double complex)*D);
		np[i] = malloc(sizeof(int)*D);//next and previous site arrays
		nm[i] = malloc(sizeof(int)*D);
		}
	initrand(seed);//initializing rand gen
	geometry(dim, L, np, nm);//calculating next and previous site, with PBC
	#ifdef DEBUG
		for(i=0; i<dim; i++){ //geometry test
			for(j=0; j<D; j++){
				if((np[nm[i][j]][j] != i) || (nm[np[i][j]][j] != i)){
					fprintf(stderr, "Geometry test failed, (%s, %d)\n", __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}
			}
	#endif
	if(confn == 10) inizialize(iflag, dim, &ds, &dg, sconfigF, field, gauge);//starting configuration, if confn == 10 loading from the final of the 10 
	else if( (confn < 10) && (confn%2 == 0)) inizialize(iflag, dim, &ds, &dg, sconfigA, field, gauge);//using in this case the last configuration which is even
	else if( (confn < 10) && (confn%2 != 0)) inizialize(iflag, dim, &ds, &dg, sconfigB, field, gauge);//using in this case the last configuration which is odd
	else inizialize(iflag, dim, &ds, &dg, sconfigF, field, gauge);//other cases, using hot or cold starting
	//energy(dim, beta, lambda, field, gauge, np);
	configA = fopen(sconfigA, "w");//saving configurations file opening, even configuration
	if(configA == NULL){
		fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	configB = fopen(sconfigB, "w");//saving configurations file opening, odd configuration (saving in two different files to prevent writing errors
	if(configB == NULL){
		fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	configF = fopen(sconfigF, "w");//saving configurations file opening, final configuration (saved only if the simulation is completed)
	if(configF == NULL){
		fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	dati = fopen(sdati, "w");//saving results file opening
	if(dati == NULL){
		fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	if((iflag == 0)||(iflag == 1)){//when not starting from a saved configuration, doing termalization steps
		for(i=0; i<iterm; i++){
			//printf("T %d\n", i);
			update(dim, ds, dg, beta, lambda, field, gauge, np, nm, acc);
			if(acc[0] < 0.3) ds *= 0.9;//autoregolation of metropolis parameters to have acceptance around 30%
			else ds *= 1.1;
			if(acc[1] < 0.3) dg *= 0.9;
			else dg *= 1.1;
			}
		}
	
	//microcanonical reversibility test, applying two times the micro step to check if it is reversibile
	#ifdef DEBUG
		for(i=0; i<dim; i++){
			for(k=0; k<N; k++) microtests[k] = field[i][k];
			micros(i, dim, beta, lambda, field, gauge, np, nm);//double microcanonical steps on scalar field
			micros(i, dim, beta, lambda, field, gauge, np, nm);
			for(k=0; k<N; k++) {
				microtests[k] -= field[i][k];//calculating the difference from the previous and new configuration, should be 0
				if(cabs(microtests[k]) > 1/pow(10, 13)) {
					fprintf(stderr, "Micro S reversibility test failed, (%s, %d)\n", __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}
			for(k=0; k<D; k++) {
				microtestg[k] = gauge[i][k];
				microg(i, k, dim, beta, lambda, field, gauge, np, nm);//double microcanonical steps on gauge field
				microg(i, k, dim, beta, lambda, field, gauge, np, nm);
				microtestg[k] -= gauge[i][k];//calculating the difference from the previous and new configuration, should be 0
				if(cabs(microtestg[k]) > 1/pow(10, 13)) {
					fprintf(stderr, "Micro G reversibility test failed, (%s, %d)\n", __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
				}
			}
	#endif
	accs = 0;
	accg = 0;
	parity = 0;
	for(i=0; i<nmeas; i++){//making nmeas measures during the markov chain
		for(j=0; j<idec; j++){//decorrelation steps before each measurement
			update(dim, ds, dg, beta, lambda, field, gauge, np, nm, acc);
			accs += acc[0];
			accg += acc[1];
			}
		temp = i%(nmeas/10);//variable used to save even, odd and final configurations
		if((i%10000) == 0) {//flush of measures every 10000 steps
			err = fflush(dati);
			if(err < 0){
				fprintf(stderr, "Error in file writing, (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			}
		measure(dim, beta, lambda, L, field, gauge, np, dati);//measurement on the configuration after microcanonical and metropolis steps
		coolfield = malloc(sizeof(double complex *)*dim);//structures allocation 
	        coolgauge = malloc(sizeof(double complex *)*dim);
        	for(j=0; j<dim; j++){
                	coolfield[j] = malloc(sizeof(double complex)*N);
                	coolgauge[j] = malloc(sizeof(double complex)*D);
                	}
		copyconf(dim, field, gauge, coolfield, coolgauge);
		updatecool(dim, beta, lambda, coolfield, coolgauge, np, nm);
		measurecool(dim, beta, lambda, L, coolfield, coolgauge, np, dati, 0);	
		updatecool(dim, beta, lambda, coolfield, coolgauge, np, nm);
                measurecool(dim, beta, lambda, L, coolfield, coolgauge, np, dati, 0); 
		updatecool(dim, beta, lambda, coolfield, coolgauge, np, nm);
                measurecool(dim, beta, lambda, L, coolfield, coolgauge, np, dati, 0); 
		updatecool(dim, beta, lambda, coolfield, coolgauge, np, nm);
                measurecool(dim, beta, lambda, L, coolfield, coolgauge, np, dati, 0); 
		updatecool(dim, beta, lambda, coolfield, coolgauge, np, nm);
                measurecool(dim, beta, lambda, L, coolfield, coolgauge, np, dati, 0); 
		updatecool(dim, beta, lambda, coolfield, coolgauge, np, nm);
                measurecool(dim, beta, lambda, L, coolfield, coolgauge, np, dati, 1);
		for(j=0; j<dim; j++) {//structures deallocation
                	free(coolfield[j]);
                	free(coolgauge[j]);
                	}
	        free(coolfield);
        	free(coolgauge);
		if((temp == 0) && (parity == 0)) {//saving the actual configuration 10 times during the total run, even configuration
			saveconfig(dim, ds, dg, field, gauge, configA);//saving configuration in even file
			log = fopen(slog, "w");//saving results file opening
			if(log == NULL){
				fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			fprintf(log, "Configuration = %d of 10\n", i/(nmeas/10)); //writing in log file the number of the saved configuration of the 10 
			fclose(log);
			parity = 1;
			}
		if((temp == 0) && (parity != 0)) {//saving the actual configuration 10 times during the total run, odd configuration
			saveconfig(dim, ds, dg, field, gauge, configB);//saving configuration in odd file
			log = fopen(slog, "w");//saving results file opening
			if(log == NULL){
				fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			fprintf(log, "Configuration = %d of 10\n", i/(nmeas/10)); //writing in log file the number of the saved configuration of the 10 
			fclose(log);
			parity = 0;
			}
		if(i == (nmeas-1)) {//saving the actual configuration 10 times during the total run, final configuration
			saveconfig(dim, ds, dg, field, gauge, configF);//saving configuration in final file
			log = fopen(slog, "w");//saving results file opening
			if(log == NULL){
				fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			fprintf(log, "Configuration = 10 of 10\n"); //writing in log file the number of the saved configuration of the 10 
			fclose(log);
			}
		}
	accs = accs/(nmeas*idec);
	accg = accg/(nmeas*idec);
	log = fopen(slog, "a");//log file to save all used parameters, after the configuration number
	if(log == NULL){
		fprintf(stderr, "Error in file opening, (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
		}
	fprintf(log, "Seed = %d\n", seed); //writing acceptance of metropolis on both fields
	fprintf(log, "N meas, term steps, decorr steps, init flag = %d, %d, %d, %d\n", nmeas, iterm, idec, iflag); //writing acceptance of metropolis on both fields
	fprintf(log, "Beta, lambda = %.2lf, %.2lf\n", beta, lambda); //writing acceptance of metropolis on both fields
	fprintf(log, "Metro S param, G param = %lf, %lf\n", ds, dg); //writing acceptance of metropolis on both fields
	for(i=0; i<D; i++) fprintf(log, "L[%d] = %d\n", i, L[i]);
	fprintf(log, "Metro S acc, G acc = %lf, %lf\n", accs, accg); //writing acceptance of metropolis on both fields
	//printf("Final metropolis parameters = S %lf, 	G %lf\n", ds, dg); //writing acceptance parameters of metropolis on both fields
	for(i=0; i<dim; i++) {//structures deallocation
		free(field[i]);
		free(gauge[i]);
		free(np[i]);
		free(nm[i]);
		}
	free(field);
	free(gauge);
	free(np);
	free(nm);
	fclose(configA);//files closing
	fclose(configB);//files closing
	fclose(configF);//files closing
	fclose(dati);
	fclose(log);
	return 0;
	}
