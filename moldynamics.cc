#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<vector>
#include<ctime>
#define random (double)rand()/(RAND_MAX +1.0)
	
	using namespace std;

struct Particle
{
	double xinit, yinit, xprev, yprev, xcur, ycur, xnew, ynew, vx, vy;
};

const int N = 4; //The number of particle in a slice
double energy(vector<Particle>, double);
double temp(vector<Particle>);
void r2(vector<Particle>,double &);
void heatup(vector<Particle> &);
void initialize(vector<Particle> &, double, double);
void initialize_lattice(vector<Particle> &, double, double);
void move(vector<Particle> &, double, double);
void relabel(vector<Particle> &);
main(){
	srand(time(NULL));
	double L = 4.0; //Size of the Box
	double dt = 0.01; //Time step
	double itmax = 100000; 	
	//Create N number of particles
	vector<Particle> particle(N*N);	
	
	//Initialize the positions and the velocities of particles
	//initialize(particle, dt, L);
	initialize_lattice(particle, dt, L);
	double vx[N*N] = {0};
	double vy[N*N] = {0};
	double v[N*N]  = {0};
	double r2mag = 0.0;
	//Evolve the system
	for(int it = 0; it < itmax; it++ ){	
		//Move the particles one step
		move(particle, dt, L);
	
		//Relabel the current position as previous
		relabel(particle);				
		
		//Printout the physical quantities
		/*int n = it%10;
		if(n==0){
		char filename[256];
		sprintf(filename, "position/%03d.dat", it/10);
		ofstream outfile(filename);
		for(int i = 0; i < N*N; i++){
			outfile << particle[i].xcur << " " << particle[i].ycur << endl;
		}
		}*/
		int m = it%100;
		if(m == 0){
			//Heat up the system
			heatup(particle);
		}
		r2(particle, r2mag);
		cout << it*dt <<"  " << energy(particle, L) << " " << temp(particle) << " " << r2mag << " " << particle[6].xcur << " "<< particle[6].ycur<<endl;
		/*if(it>itmax*1/5){
			for(int i = 0; i < N*N; i++){
				 vx[i] += particle[i].vx;
				 vy[i] += particle[i].vy;
				 v[i]  += sqrt(particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy );
			}
		}*/
	}
//	ofstream outfile("velocity.dat");
//	for(int i = 0; i < N*N; i++){
//		outfile << vx[i]/(itmax/4*5)<< "  " << v[i]/(itmax/4*5) << endl;
//	}
}

void initialize(vector<Particle> &particle, double dt, double L){
	double v0 = 1.0;  //Initial velocity amplitude	
	double del = 0.5;
	int m, n;
	for(int i = 0; i < N*N; i++){
			m = int(i/N);
			n = int(i%N);
			particle[i].xcur = 2.0*(m+1) + 2*(random-0.5)*del;
			particle[i].ycur = 2.0*(n+1) + 2*(random-0.5)*del;
			particle[i].xinit = particle[i].xcur;
			particle[i].yinit = particle[i].ycur;
			particle[i].vx = 2*(random-0.5)*v0; 
			particle[i].vy = 2*(random-0.5)*v0; 
			particle[i].xprev = particle[i].xcur - particle[i].vx*dt;
			particle[i].yprev = particle[i].ycur - particle[i].vy*dt;
	}
}
void initialize_lattice(vector<Particle> &particle, double dt, double L){
	double v0 = 0.0001;  //Initial velocity amplitude	
	int m, n;
	for(int i = 0; i < N*N; i++){
			m = int(i/N);
			n = int(i%N);
			particle[i].xcur = m;
			particle[i].ycur = n;
			particle[i].xinit = particle[i].xcur;
			particle[i].yinit = particle[i].ycur;
			particle[i].vx = pow(-1,rand()%2)*v0; 
			particle[i].vy = pow(-1,rand()%2)*v0; 
			particle[i].xprev = particle[i].xcur - particle[i].vx*dt;
			particle[i].yprev = particle[i].ycur - particle[i].vy*dt;
	}

}
void move(vector<Particle> &particle, double dt, double L){
	double fx[N*N];
	double fy[N*N];
	for(int k = 0; k < N*N; k++){
		for(int i = 0; i < N*N; i++){
			double r, rx, ry;
			double Fx = 0.0;
			double Fy = 0.0;
			for(int j = 0; j < N*N; j++){
				if(j==i) continue;	
				rx = particle[i].xcur - particle[j].xcur;  	
				ry = particle[i].ycur - particle[j].ycur;
				if(abs(rx) > 0.5*L) rx *= -L/abs(rx) + 1.0;
				if(abs(ry) > 0.5*L) ry *= -L/abs(ry) + 1.0;
				r = sqrt(rx*rx + ry*ry);
				if(r > 3.0) continue;
				//if(r < 0.1) continue;
				Fx += 24*(2.0/pow(r,13)-1.0/pow(r,7))*rx/r;
				Fy += 24*(2.0/pow(r,13)-1.0/pow(r,7))*ry/r;
			}
			//cout << i << " " << Fx <<" "<< Fy << endl;
			fx[i] = Fx;
			fy[i] = Fy;
		}
		particle[k].xnew = (2*particle[k].xcur - particle[k].xprev + fx[k]*dt*dt);
		if(particle[k].xnew > L) particle[k].xnew -= L; 	
		if(particle[k].xnew < 0) particle[k].xnew += L; 	
		particle[k].ynew = (2*particle[k].ycur - particle[k].yprev + fy[k]*dt*dt); 	
		if(particle[k].ynew > L) particle[k].ynew -= L; 	
		if(particle[k].ynew < 0) particle[k].ynew += L; 	
		particle[k].vx = (particle[k].xnew - particle[k].xprev)/(2*dt);	
		particle[k].vy = (particle[k].ynew - particle[k].yprev)/(2*dt);	
	
		//cout << k <<" " << particle[k].xnew << " "<< particle[k].xprev << " " << particle[k].vx << endl;
	}
}
void relabel(vector<Particle> &particle){
	for(int i = 0; i < N*N; i++){
		particle[i].xprev = particle[i].xcur; 
		particle[i].yprev = particle[i].ycur; 
		particle[i].xcur = particle[i].xnew; 
		particle[i].ycur = particle[i].ynew; 
	}
}
double energy(vector<Particle> particle, double L){
	double KE = 0.0;
	double PE = 0.0;
	double rx, ry, r; 
	for(int i = 0; i < N*N; i++){
		for(int j = 0; j < N*N; j++){
			if(i==j) continue;
			rx = particle[i].xcur - particle[j].xcur;  	
			ry = particle[i].ycur - particle[j].ycur;
			if(abs(rx) > 0.5*L) rx *= -L/abs(rx) + 1.0;
			if(abs(ry) > 0.5*L) ry *= -L/abs(ry) + 1.0;
		  	r = sqrt(rx*rx + ry*ry);
			if(r > 3.0) continue;
			//if(r < 0.1) continue;
			PE += 2*(1.0/pow(r,12)-1.0/pow(r,6));	
		}
		KE += 0.5*(particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy);
	}
	return KE + PE;
}
double temp(vector<Particle> particle){
	double KE = 0.0;
	for(int i = 0; i < N*N; i++){
		KE += 0.5*(particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy); 
	}	
	return KE/(N*N-1);	
}

void heatup(vector<Particle> &particle){
	for(int i = 0; i < N*N; i++){
		particle[i].vx *= 2;
		particle[i].vy *= 2;
	}	
}

void r2(vector<Particle> particle, double &r2mag){
	double delx, dely;
	int i = 5;
	//for(int i = 0; i < N*N; i++){
		delx = (particle[i].xcur - particle[i].xinit);		
		dely = (particle[i].ycur - particle[i].yinit);		
		r2mag = (delx*delx + dely*dely); 
	//}
}


