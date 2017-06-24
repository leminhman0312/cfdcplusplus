//Solving the 2D linear convection equation 

#include <iostream>
#include <cmath>
#include "engine.h"
#include <string.h>
using namespace std;

#define BUFSIZE 256


int main(){
	
	//Declaring variables
	int tmax = 2.0;
	int xmax = 2.0;
	int xmin = 0.0;
	int ymax = 2.0;
	int ymin = 0.0;
	int nx = 80; //x steps
	int ny = 80; // y steps
	int nt = 100; //time steps 	
	int c = 1; //constant
	double dx = (xmax-xmin)/double(nx-1);
	double dy = (ymax-ymin)/double(ny-1);
	double sigma = 0.2;
	double dt = sigma*dx;
	double x[nx], y[ny]; //spatial vectors X and Y
	double u[ny][nx]; //solution matrix
	double u0[ny][nx]; //initial matrix
	


	//Setting up initial conditions

	for (int j = 0; j<ny; j++){
		y[j] = j*dy;
		for (int i = 0; i <nx;i++){
			x[i] = i*dx;
			if((x[i]>= 0.5 && x[i]<1.0) && (y[j] >=0.5 && y[j] < 1.0)){
				u0[j][i] = 2.0;
			}
			else{
				u0[j][i] = 1.0;
			}
		}
	}

	

	//Plot the initial conditions (it works!)

	Engine *ep;

	//Meshgrid 
	mxArray *x_pointer = NULL;
	mxArray *y_pointer = NULL;
	mxArray *u_pointer = NULL;
	mxArray *nt_pointer = NULL; //for looping in MATLAB later
	mxArray *u0_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//Test if Matlab is opened or not
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return EXIT_FAILURE;
	}

	//X VECTOR
	x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
	memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
	engPutVariable(ep, "X_matlab", x_pointer);

	y_pointer = mxCreateDoubleMatrix(1,ny,mxREAL);
	memcpy((void *)mxGetPr(y_pointer), (void *)y, sizeof(y));
	engPutVariable(ep, "Y_matlab", y_pointer);

	u0_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
	memcpy((void *)mxGetPr(u0_pointer), (void *)u0, sizeof(u0));
	engPutVariable(ep, "U0_matlab", u0_pointer);




	//MESH GRID

	engEvalString(ep,"[X Y] = meshgrid(X_matlab,Y_matlab);");
	

	//CALCULATIONS FOR THE ACTUAL SOLUTION

	for (int n = 0;n<=nt;n++){
		//reassign to a holder U0
		for (int j = 1; j<=ny-1;j++){
			for (int i = 1; i<=nx-1;i++){
				u[j][i] = u0[j][i]; //copy
			}
		}

		//Now for actual calculations
		for (int k = 1; k<=ny-1;k++){
			for (int l = 1; l<=nx-1;l++){
				u0[k][l] = u[k][l] - (c*dt/dx*(u[k][l] - u[k][l-1]))-(c*dt/dy*(u[k][l]-u[k-1][l]));
			}
		}

		//Now for the boundary conditions
		
		for (int w = 0; w<nx;w++){
			u[0][w] = 1; //TOP 
			u[w][ny-1] = 1; //RIGHT
		}

		for (int s = 0;s<ny;s++){
			u[s][0] = 1; //LEFT
			u[nx-1][s] =1; //BOTTOM
		}

		
		
		//FOR ANIMATED PLOT, uncomment this
		/*
		u_pointer = mxCreateDoubleMatrix(ny,nx,mxREAL);
		memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
		engPutVariable(ep, "U_matlab", u_pointer);
		engEvalString(ep,"fig1= surf(X_matlab,Y_matlab,U_matlab);");
		engEvalString(ep,"pause(0.1);");
		*/ 
				
	}

	
	u_pointer = mxCreateDoubleMatrix(ny,nx,mxREAL);
	memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
	engPutVariable(ep, "U_matlab", u_pointer);
	engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
	engEvalString(ep,"fig1 = surf(X_matlab,Y_matlab,U_matlab);");
	engEvalString(ep,"view(0,90);");
	engEvalString(ep,"colorbar('southoutside');");
	engEvalString(ep,"xlabel('X VALUES');");
	engEvalString(ep,"ylabel('Y VALUES');");
	engEvalString(ep,"zlabel('U VELOCITY');");
	engEvalString(ep,"title('2D LINEAR CONVECTION');");
	engEvalString(ep,"saveas(fig1, '2D_LINEAR_CONVECTION.png');");

	
	
	

	return 0;
}