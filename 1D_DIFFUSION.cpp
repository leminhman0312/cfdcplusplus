#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine.h"
#include <cmath>
using namespace std;


#define BUFSIZE 256


//SOLVING THE 1D LINEAR DIFFUSION EQUATION

// dU/dt = Nu * d^2(U)/dx^2


int main(){

	
	///////////////////////////////////////////////////
	////////////////Start of C++ /////////////////////

	//Declare basic variables


	int xmax = 2;
	int xmin = 0;
	int nx = 41;
	double dx = double((xmax-xmin)) /(nx-1);


	int nt = 20;
	double nu = 0.3; //viscosity
	double sigma = 0.2; // a parameter
	double dt = pow(dx,2.0)*sigma/nu;

	//cout << dt << endl; 
	
	
	double x[nx]; //space vectors
	double u0[nx]; //U at time = 0
	double u[nx]; //Real U[time][space], u[i][j]


		

	//Initial Condition.  Step function

	for (int j = 0; j < nx; j++){
		x[j] = j*dx;
		if (x[j]<0.5 || x[j] > 1.0){
			u0[j] = 1.0;
		}
		else{
			u0[j] = 2.0;
		}

	}


	/*
	for (int n = 0; n<=nt; n++){
		//assigning value at previous time
		//using u0 as a holder
		for (int j = 0; j<nx;j++){
			u[j] = u0[j];
		}
		for (int i =1; i< nx-1; i++){
			u0[i]=u[i]+nu*dt/(pow(dx,2.0))*(u[i+1]-2*u[i]+u[i-1]);

		}

	}

	//////////////  END C++ //////////////////////////////////////////

	//////////////////////////////////////////////////////////////////
	//////////////START OF MATLAB PROGRAMMING////////////////////////

	
	
	//Declare basic pointers
	Engine *ep;
	mxArray *u_pointer = NULL;
	mxArray *x_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//Putting C array into MATLAB array
	
	//X VECTOR
	x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
	memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
	engPutVariable(ep, "X_matlab", x_pointer);

	//U VECTOR

	u_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
	memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
	engPutVariable(ep, "U_matlab", u_pointer);


	//PLOTTING U
	
	
	
	engEvalString(ep,"k = plot(X_matlab,U_matlab);");
	engEvalString(ep,"grid on");
	engEvalString(ep,"box on");
	engEvalString(ep,"title('1D DIFFUSION')");
	engEvalString(ep,"xlabel('Space')");
	engEvalString(ep,"ylabel('Velocity')");
	engEvalString(ep,"saveas(k,'1D_DIFFUSION.png')");

	
	mxDestroyArray(u_pointer);
	mxDestroyArray(x_pointer);
	engEvalString(ep, "close;"); 
	
	
	return EXIT_SUCCESS;
	*/ 





	//FOR VIDEOS

	//Declare basic pointers
	Engine *ep = engOpen(NULL);
	mxArray *u_pointer = NULL;
	mxArray *x_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//OPEN VIDEOS

	engEvalString(ep,"v=VideoWriter('1D_DIFFUSION.avi');");
	engEvalString(ep,"open(v)");

	//LOOP WITH VIDEOS

	for (int n = 0; n<=nt; n++){
		//assigning value at previous time
		//using u0 as a holder
		for (int j = 0; j<nx;j++){
			u[j] = u0[j];
		}
		for (int i =1; i< nx-1; i++){
			u0[i]=u[i]+nu*dt/(pow(dx,2.0))*(u[i+1]-2*u[i]+u[i-1]);

		}
		//X VECTOR
		x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
		memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
		engPutVariable(ep, "X_matlab", x_pointer);

		//U VECTOR

		u_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
		memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
		engPutVariable(ep, "U_matlab", u_pointer);

		//ANIMATING

		engEvalString(ep,"plot(X_matlab,U_matlab);");
		engEvalString(ep,"grid on");
		engEvalString(ep,"box on");
		engEvalString(ep,"title('1D DIFFUSION')");
		engEvalString(ep,"xlabel('Space')");
		engEvalString(ep,"ylabel('Velocity')");
		engEvalString(ep,"M=getframe(gcf);");
		engEvalString(ep,"writeVideo(v,M)");

	}
	engEvalString(ep,"close(v)");

}


















