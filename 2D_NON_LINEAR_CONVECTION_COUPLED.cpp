//Solving the 2D NON LINEAR CONVECTION EQUATIOn

// Dv/Dt + u*(Dv/Dx) + v(Dv/Dy) = 0 

#include <iostream>
#include "engine.h"
#include <string.h>
#include <cmath>

using namespace std;

#define BUFSIZE 256

int main(){

	//Declare Variables

	int nx = 101;
	int ny = 101;
	int nt = 80;
	int c = 1;
	double sigma = 0.2;
	double dx = 2/double((nx-1));
	double dy = 2/double((ny-1));
	double dt = sigma*dx;
	double x[nx], y[ny]; //spatial vectors X and Y
	

	double u[ny][nx], v[ny][nx]; //solution matrix, 2 of them, 1 for X. 1 for Y

	double u0[ny][nx], v0[ny][nx]; //initial matrix

	//Setting up boundary conditions (all = 1)
	//SETUP BOUNDARY CONDITION, USING LOOPS


	for (int i = 0; i<nx;i++){
		for (int j = 0; j<ny;j++){
			u[j][i] =1;
			u0[j][i]=1;
			v[j][i] =1;
			v0[j][i]=1;
		}
	}


	//Setting up Initial Conditions

	//For U
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

	//For V
	for (int j = 0; j<ny; j++){
		y[j] = j*dy;
		for (int i = 0; i <nx;i++){
			x[i] = i*dx;
			if((x[i]>= 0.5 && x[i]<1.0) && (y[j] >=0.5 && y[j] < 1.0)){
				v0[j][i] = 2.0;
			}
			else{
				v0[j][i] = 1.0;
			}
		}
	}



	
	
	
	
	//for testing IC
	//engEvalString(ep,"fig1 = surf(X_matlab,Y_matlab,U0_matlab);");
	//engEvalString(ep,"saveas(fig1,'2D_NON_LINEAR_IC_U.png');");
	
	/*

	//CALCULATIONS

	for (int n = 0; n<=nt; n++){
		//Copy Loops into U0, V0;
		//Reassign values
		
		//COPY FOR U
		for (int j = 1; j<=ny-1;j++){
			for (int i = 1; i<=nx-1;i++){
				u[j][i] = u0[j][i]; //copy
			}
		}

		//COPY FOR V
		for (int k = 1; k<=ny;k++){
			for (int l = 1; l<=nx;l++){
				v[k][l] = v0[k][l]; //copy
			}
		}

		//MAIN LOOP
		// LOOP FOR real U and V
		for (int yiter = 1; yiter<=ny-1;yiter++){
			for (int xiter = 1; xiter<=nx-1;xiter++){
				u0[yiter][xiter] = u[yiter][xiter]-(u[yiter][xiter]*(c*dt/dx)*(u[yiter][xiter]-u[yiter][xiter-1]))-(v[yiter][xiter]*(c*dt/dy)*(u[yiter][xiter]-u[yiter-1][xiter]));
				v0[yiter][xiter] = v[yiter][xiter]-(u[yiter][xiter]*(c*dt/dx)*(v[yiter][xiter]-v[yiter][xiter-1]))-(v[yiter][xiter]*(c*dt/dy)*(v[yiter][xiter]-v[yiter-1][xiter]));
			}
		}	

	}


	
	//PLUG U AND V INTO MATLAB

	Engine *ep;
	mxArray *x_pointer = NULL;
	mxArray *y_pointer = NULL;
	mxArray *u_pointer = NULL;
	mxArray *v_pointer = NULL;
	mxArray *u0_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];


	//Test if Matlab is opened or not
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return EXIT_FAILURE;
	}

	x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
	memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
	engPutVariable(ep, "X_matlab", x_pointer);

	y_pointer = mxCreateDoubleMatrix(1,ny,mxREAL);
	memcpy((void *)mxGetPr(y_pointer), (void *)y, sizeof(y));
	engPutVariable(ep, "Y_matlab", y_pointer);

	u0_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
	memcpy((void *)mxGetPr(u0_pointer), (void *)u0, sizeof(u0));
	engPutVariable(ep, "U0_matlab", u0_pointer);

	
	u_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
	memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
	engPutVariable(ep, "U_matlab", u_pointer);


	v_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
	memcpy((void *)mxGetPr(v_pointer), (void *)v, sizeof(v));
	engPutVariable(ep, "V_matlab", v_pointer);

	//MESH GRID
	engEvalString(ep,"[X Y] = meshgrid(X_matlab,Y_matlab);");


	//PLOTTING

	engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
	engEvalString(ep,"fig2 = surf(X_matlab,Y_matlab,U_matlab);");
	engEvalString(ep,"colormap(jet);");
	engEvalString(ep,"shading interp;");
	engEvalString(ep,"view(0,90);");
	engEvalString(ep,"colorbar('southoutside');");
	engEvalString(ep,"xlabel('X VALUES');");
	engEvalString(ep,"ylabel('Y VALUES');");
	engEvalString(ep,"zlabel('U VELOCITY');");
	engEvalString(ep,"title('2D NON LINEAR CONVECTION U ');");
	engEvalString(ep,"saveas(fig2, '2D_NON_LINEAR_CONVECTION_X.png');");


	engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
	engEvalString(ep,"fig3 = surf(X_matlab,Y_matlab,V_matlab);");
	engEvalString(ep,"colormap(jet);");
	engEvalString(ep,"shading interp;");
	engEvalString(ep,"view(0,90);");
	engEvalString(ep,"colorbar('southoutside');");
	engEvalString(ep,"xlabel('X VALUES');");
	engEvalString(ep,"ylabel('Y VALUES');");
	engEvalString(ep,"zlabel('V VELOCITY');");
	engEvalString(ep,"title('2D NON LINEAR CONVECTION V');");
	engEvalString(ep,"saveas(fig3, '2D_NON_LINEAR_CONVECTION_Y.png');");

	
	*/ 

	// WITH VIDEOS

	// MATLAB variables

	Engine *ep = engOpen(NULL);
	mxArray *x_pointer = NULL;
	mxArray *y_pointer = NULL;
	mxArray *u_pointer = NULL;
	mxArray *v_pointer = NULL;
	mxArray *u0_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//FOR VIDEO

	//engEvalString(ep,"vx = VideoWriter('2D_NON_LINEAR_CONVECTION_X.avi');");
	//engEvalString(ep,"open(vx)");

	engEvalString(ep,"vy = VideoWriter('2D_NON_LINEAR_CONVECTION_Y.avi');");
	engEvalString(ep,"open(vy)");
	
	//CALCULATIONS with videos

	for (int n = 0; n<=nt; n++){
		//Copy Loops into U0, V0;
		//Reassign values
		
		//COPY FOR U
		for (int j = 1; j<=ny-1;j++){
			for (int i = 1; i<=nx-1;i++){
				u[j][i] = u0[j][i]; //copy
			}
		}

		//COPY FOR V
		for (int k = 1; k<=ny;k++){
			for (int l = 1; l<=nx;l++){
				v[k][l] = v0[k][l]; //copy
			}
		}

		//MAIN LOOP
		// LOOP FOR real U and V
		for (int yiter = 1; yiter<=ny-1;yiter++){
			for (int xiter = 1; xiter<=nx-1;xiter++){
				u0[yiter][xiter] = u[yiter][xiter]-(u[yiter][xiter]*(c*dt/dx)*(u[yiter][xiter]-u[yiter][xiter-1]))-(v[yiter][xiter]*(c*dt/dy)*(u[yiter][xiter]-u[yiter-1][xiter]));
				v0[yiter][xiter] = v[yiter][xiter]-(u[yiter][xiter]*(c*dt/dx)*(v[yiter][xiter]-v[yiter][xiter-1]))-(v[yiter][xiter]*(c*dt/dy)*(v[yiter][xiter]-v[yiter-1][xiter]));
			}
		}

		x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
		memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
		engPutVariable(ep, "X_matlab", x_pointer);

		y_pointer = mxCreateDoubleMatrix(1,ny,mxREAL);
		memcpy((void *)mxGetPr(y_pointer), (void *)y, sizeof(y));
		engPutVariable(ep, "Y_matlab", y_pointer);

		u_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
		memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
		engPutVariable(ep, "U_matlab", u_pointer);

		v_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
		memcpy((void *)mxGetPr(v_pointer), (void *)v, sizeof(v));
		engPutVariable(ep, "V_matlab", v_pointer);


		/*
		//PLOTTING X

		engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
		engEvalString(ep,"fig2 = surf(X_matlab,Y_matlab,U_matlab);");
		engEvalString(ep,"colormap(jet);");
		engEvalString(ep,"shading interp;");
		engEvalString(ep,"view(30,60);");
		engEvalString(ep,"colorbar('southoutside');");
		engEvalString(ep,"xlabel('X VALUES');");
		engEvalString(ep,"ylabel('Y VALUES');");
		engEvalString(ep,"zlabel('U VELOCITY');");
		engEvalString(ep,"title('2D NON LINEAR CONVECTION U ');");
		engEvalString(ep,"M = getframe(gcf);");
		engEvalString(ep,"writeVideo(vx,M)");
		*/ 


		//PLOTTING Y 

		engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
		engEvalString(ep,"fig3 = surf(X_matlab,Y_matlab,V_matlab);");
		engEvalString(ep,"colormap(jet);");
		engEvalString(ep,"shading interp;");
		engEvalString(ep,"view(30,60);");
		engEvalString(ep,"colorbar('southoutside');");
		engEvalString(ep,"xlabel('X VALUES');");
		engEvalString(ep,"ylabel('Y VALUES');");
		engEvalString(ep,"zlabel('V VELOCITY');");
		engEvalString(ep,"title('2D NON LINEAR CONVECTION V');");
		engEvalString(ep,"M = getframe(gcf);");
		engEvalString(ep,"writeVideo(vy,M)");


		
		}

	engEvalString(ep,"close(vy)");

	return 0;
}