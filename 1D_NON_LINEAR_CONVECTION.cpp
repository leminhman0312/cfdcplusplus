#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine.h"

using namespace std;


#define BUFSIZE 256


//SOLVING THE 1D LINEAR CONVECTION

// du/dt + U*du/dx = 0


int main(){

	
	///////////////////////////////////////////////////
	////////////////Start of C++ /////////////////////

	//Declare basic variables
	int tmax = 24;
	int tmin = 0;
	int nt = 25;
	//float dt = float((tmax-tmin))/float((nt-1));
	float dt = 0.025;

	int xmax = 2;
	int xmin = 0;
	int nx = 41;
	float dx = float((xmax-xmin)) /(float((nx-1)));
	int c = 1;
	
	double x[nx]; //space vectors
	double u0[nx]; //U at time = 0
	double u[nx]; //u[time][space], u[i][j]

	float alpha = c*dt/dx;

		

	//Initial Condition.  Step function

	for (int j = 0; j < nx; j++){
		x[j] = j*dx;
		if (x[j]>=0.5 && x[j] <=1.0){
			u0[j] = 2.0;
		}
		else{
			u0[j] = 1.0;
		}

	}

	/* NO VIDEOS

	//CALCULATIONS

	for (int i = 0; i<=nt; i++){
		for (int j = 0; j<nx;j++){
			u[j] = u0[j];
		}
		for (int j =1; j< nx; j++){
			u0[j] = u[j] - (u[j]*dt/dx)*(u[j]-u[j-1]);
		}

	}

	for (int i = 0; i< nx; i++){
		//cout << x[i] << "\t\t\t" << u[i] << endl;

	}





	//////////////  END C++ //////////////////////////////////////////

	//////////////////////////////////////////////////////////////////
	//////////////START OF MATLAB PROGRAMMING////////////////////////


	
	//Declare basic pointers
	Engine *ep;
	mxArray *u_pointer = NULL;
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


	//PLOTTING U
	
	
	engEvalString(ep,"g = plot(X_matlab,U_matlab);");
	engEvalString(ep,"grid on");
	engEvalString(ep,"box on");
	engEvalString(ep,"title('1D LINEAR CONVECTION')");
	engEvalString(ep,"xlabel('Space')");
	engEvalString(ep,"ylabel('Velocity')");
	engEvalString(ep,"saveas(g,'1D_NON_LINEAR_CONVECTION.png')");

	
	mxDestroyArray(u_pointer);
	mxDestroyArray(x_pointer);
	engEvalString(ep, "close;");

	
	
	return EXIT_SUCCESS;
	
	*/ 


	//VIDEO PARTS

	//Declare basic pointers
	Engine *ep = engOpen(NULL);
	mxArray *u_pointer = NULL;
	mxArray *x_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//OPEN VIDEOS

	engEvalString(ep,"v = VideoWriter('1D_NON_LINEAR_CONVECTION.avi');");
	engEvalString(ep,"open(v)");

	//CALCULATIONS with VIDEO

	for (int i = 0; i<=nt; i++){
		for (int j = 0; j<nx;j++){
			u[j] = u0[j];
		}
		for (int j =1; j< nx; j++){
			u0[j] = u[j] - (u[j]*dt/dx)*(u[j]-u[j-1]);
		}

		//X VECTOR
		x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
		memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
		engPutVariable(ep, "X_matlab", x_pointer);
		//U VECTOR
		u_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
		memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
		engPutVariable(ep, "U_matlab", u_pointer);
		//PLOTTING U
		engEvalString(ep,"plot(X_matlab,U_matlab);");
		engEvalString(ep,"grid on");
		engEvalString(ep,"box on");
		engEvalString(ep,"title('1D LINEAR CONVECTION')");
		engEvalString(ep,"xlabel('Space')");
		engEvalString(ep,"ylabel('Velocity')");
		engEvalString(ep,"M = getframe(gcf);");
		engEvalString(ep,"writeVideo(v,M)");
	}

	engEvalString(ep,"close(v)");

	return 0;
}


















