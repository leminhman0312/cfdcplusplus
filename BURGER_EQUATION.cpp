//Solving Burgers Equation 
// dU/dt + u(Du/Dx) = v(D^2(u)/Dx^2)
// Combination of Non Linear Convection and Diffusion

#include <iostream>
#include <cmath>
#include "engine.h"
#include <stdio.h>
#include <string.h>
using namespace std;

#define BUFSIZE 256


int main(){


	/////////////////////////////////////////////////////////////////////

	//START OF c++ ///////////////////////////////////////////////////

	
	//Declare variables

	int nx = 101; //nodes number
	int nt = 100; //time steps
	double nu = 0.07; //viscosity 
	const double pi = 4*atan(1);
	double dx = 2*pi/(nx-1);
	double dt = dx*nu;

	//Setting up initial conditions

	double x[nx], phi[nx], dphi[nx]; //numerical variables
	double phir[nx], dphir[nx]; //analytical variables (real)
	double u[nx], u0[nx], u00[nx]; //velocity, velocity at t = 0;
	double ur[nx]; // analytical solution

	//Calculate ICs

	for (int i = 0; i<nx; i++){
		x[i] = i*dx;
		phi[i] = exp(-1.*pow(x[i],2.)/4./nu)+exp(-1.*pow((x[i]-2*pi),2.)/4./nu);
		dphi[i]=-0.5/nu*x[i]*exp(-1.*pow(x[i],2.)/4./nu)-0.5/nu*(x[i]-2*pi)*exp(-1.*pow((x[i]-2*pi),2.)/4./nu);
		u00[i]=u0[i]=-2.0*nu*dphi[i]/phi[i]+4.0; //initial U
	}

	//Boundary condition i-1 and i+1 terms

	// u(0) = u(2pi).  So as you go the the right, U wraps around to the front of frame

	// So say you start at 0 to nx.  1 to nx-1 is in the loop => OK
	// point 0 and point (nx) are the boundary conditions. 
	// This loop is to make sure that when we reach those points

	
	int ip1[nx], im1[nx];

	for (int i = 0; i<nx;i++){
		ip1[i] = i+1;
		im1[i] = i-1;	
	}
	ip1[nx-1] = 0;

	im1[0] = nx-1;
	




	//Populate anylytical Solution 

	double tid=(nt-1)*dt;
	for (int i = 0; i <= nx; i++) {

		phir[i]=exp(-1.*pow((x[i]-4.*tid),2.)/4./nu/(tid+1.))+exp(-1.*pow((x[i]-4.*tid-2*pi),2.)/4./nu/(tid+1.));

		dphir[i]=-0.5/nu/(tid+1.)*(x[i]-4.*tid)*exp(-1.*pow((x[i]-4.*tid),2.)/4./nu/(tid+1.))-0.5/nu/(tid+1.)*(x[i]-4*tid-2*pi)*exp(-1.*pow((x[i]-4.*tid-2*pi),2.)/4./nu/(tid+1.));

		ur[i]=-2.0*nu*dphir[i]/phir[i]+4.0;

	}

	//Calculate numerical solution 

	for (int n = 0; n<=nt;n++){
		for (int j = 0; j<nx;j++){
			u[j] = u0[j];
		}
		
		for (int i = 0; i<nx;i++){
			u0[i]=u[i]-u[i]*dt/dx*(u[i]-u[im1[i]])+nu*dt/(pow(dx,2.0))*(u[ip1[i]]-2*u[i]+u[im1[i]]);
		}			
	}


	//Test printing out results

	for (int i = 0; i < nx; i++){
		cout << ur[i] << "\t" << endl;
	}


	/////////////////////////////////////////////////////////////////////

	//START OF MATLAB ///////////////////////////////////////////////////


	//Declare pointers

	Engine *ep;
	mxArray *u_pointer = NULL;
	mxArray *u_real_pointer = NULL;
	mxArray *x_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//Test if Matlab is opened or not
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return EXIT_FAILURE;
	}

	//Putting C array into MATLAB array

	//X VECTOR
	x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
	memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
	engPutVariable(ep, "X_matlab", x_pointer);

	//U VECTOR

	u_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
	memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
	engPutVariable(ep, "U_matlab", u_pointer);

	//U REAL VECTOR

	u_real_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
	memcpy((void *)mxGetPr(u_real_pointer), (void *)ur, sizeof(ur));
	engPutVariable(ep, "U_real_matlab", u_real_pointer);


	//PLOTTING//
	engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
	engEvalString(ep,"b = plot(X_matlab,U_matlab,'xr');");
	engEvalString(ep,"axis([0 2*pi 0 10]);");
	engEvalString(ep,"grid on");
	engEvalString(ep,"box on");
	engEvalString(ep,"hold on");
	engEvalString(ep,"b = plot(X_matlab,U_real_matlab,'-b');");
	engEvalString(ep,"title('Burgers Equation')");
	engEvalString(ep,"xlabel('Space')");
	engEvalString(ep,"ylabel('Velocity')");
	engEvalString(ep,"legend('Numerical','Analytical');");
	engEvalString(ep,"saveas(b,'BURGERS.png')");

	
	mxDestroyArray(u_pointer);
	mxDestroyArray(x_pointer);
	engEvalString(ep, "close;"); 


	return 0;
}
