//SOLVING THE 2D DIFFUSION EQUATION
// DU/DT = v*(D^2(U)/(DX^2)) + v*(D^U/DY^2)


#include <stdio.h>
#include <iostream>
#include <string.h>
#include "engine.h"
#include <cmath>
using namespace std;
#define BUFSIZE 256


double centralDiff(double front, double middle, double end){
	return front-2*middle+end;

}



int main(){
	//SET UP VARIABLES
	int nx = 10;
	int ny = 10;
	int nt = 100;
	double nu = 0.05;
	double dx = 2/double((nx-1));
	double dy = 2/double((ny-1));
	double sigma = 0.25;
	double dt = (sigma*dx*dy)/(nu);
	double x[nx], y[ny]; //spatial vectors X and Y
	double u[ny][nx]; //Real solution at n+1
	double u0[ny][nx]; //Initial U, storage
	double constX = (nu*dt)/(pow(dx,2));
	double constY = (nu*dt)/(pow(dy,2));


	//SETUP BOUNDARY CONDITION, USING LOOPS


	for (int i = 0; i<nx;i++){
		for (int j = 0; j<ny;j++){
			u[j][i] =1;
			u0[j][i]=1;
		}
	}




	//SET UP INITIAL CONDITIONS 

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


	/*NO VIDEOS
	//CALCULATIONS

	for (int n = 0; n<nt; n++){
		//Reassign values		
		//COPY FOR U

		for (int j = 1; j<ny-1;j++){
			for (int i = 1; i<nx-1;i++){
				u[j][i] = u0[j][i]; //copy
			}
		}

		//Main loop, dont count first and last elements
		for (int yiter = 1; yiter<ny-1;yiter++){
			for (int xiter = 1; xiter<nx-1;xiter++){
				u0[yiter][xiter] = u[yiter][xiter] + (constX*(centralDiff(u[yiter][xiter+1],u[yiter][xiter],u[yiter][xiter-1])))+(constY*(centralDiff(u[yiter+1][xiter],u[yiter][xiter],u[yiter-1][xiter])));
			}
		}
	}

	for (int i = 0; i<nx;i++){
		for (int j = 0; j<ny;j++){
			cout << u[j][i] << "\t";
		}
		cout <<"\n";
	}
	

	//PLOTTING

	
	
	Engine *ep;
	mxArray *x_pointer = NULL;
	mxArray *y_pointer = NULL;
	mxArray *u_pointer = NULL;
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

	u_pointer = mxCreateDoubleMatrix(ny,nx,mxREAL);
	memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
	engPutVariable(ep, "U_matlab", u_pointer);

	engEvalString(ep,"[X Y] = meshgrid(X_matlab,Y_matlab);");
	

	//engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
	
	
	engEvalString(ep,"fig5 = surf(X_matlab,Y_matlab,U_matlab);");
	engEvalString(ep,"shading interp;");
	engEvalString(ep,"colormap(jet);");
	engEvalString(ep,"colorbar('southoutside');");
	engEvalString(ep,"xlabel('X VALUES');");
	engEvalString(ep,"ylabel('Y VALUES');");
	engEvalString(ep,"zlabel('U VELOCITY');");
	engEvalString(ep,"title('2D LINEAR DIFFUSION 3D');");
	engEvalString(ep,"saveas(fig5, '2D_DIFFUSION_ALL.png');");
	fgetc(stdin);

	*/ 


	//VIDEO PARTS

	Engine *ep = engOpen(NULL);
	mxArray *x_pointer = NULL;
	mxArray *y_pointer = NULL;
	mxArray *u_pointer = NULL;
	mxArray *u0_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//FOR VIDEO

	engEvalString(ep,"v = VideoWriter('2D_DIFFUSION_ALL.avi');");
	engEvalString(ep,"open(v)");

	//CALCULATIONS with VIDEO

	for (int n = 0; n<nt; n++){
		//Reassign values		
		//COPY FOR U

		for (int j = 1; j<ny-1;j++){
			for (int i = 1; i<nx-1;i++){
				u[j][i] = u0[j][i]; //copy
			}
		}

		//Main loop, dont count first and last elements
		for (int yiter = 1; yiter<ny-1;yiter++){
			for (int xiter = 1; xiter<nx-1;xiter++){
				u0[yiter][xiter] = u[yiter][xiter] + (constX*(centralDiff(u[yiter][xiter+1],u[yiter][xiter],u[yiter][xiter-1])))+(constY*(centralDiff(u[yiter+1][xiter],u[yiter][xiter],u[yiter-1][xiter])));
			}
		}

		x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
		memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
		engPutVariable(ep, "X_matlab", x_pointer);

		y_pointer = mxCreateDoubleMatrix(1,ny,mxREAL);
		memcpy((void *)mxGetPr(y_pointer), (void *)y, sizeof(y));
		engPutVariable(ep, "Y_matlab", y_pointer);

		u_pointer = mxCreateDoubleMatrix(ny,nx,mxREAL);
		memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
		engPutVariable(ep, "U_matlab", u_pointer);

		//PLOTTING

		engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure	
		engEvalString(ep,"fig5 = surf(X_matlab,Y_matlab,U_matlab);");
		engEvalString(ep,"shading interp;");
		engEvalString(ep,"colormap(jet);");
		engEvalString(ep,"colorbar('southoutside');");
		engEvalString(ep,"xlabel('X VALUES');");
		engEvalString(ep,"ylabel('Y VALUES');");
		engEvalString(ep,"zlabel('U VELOCITY');");
		engEvalString(ep,"title('2D LINEAR DIFFUSION 3D');");
		engEvalString(ep,"M = getframe(gcf);");
		engEvalString(ep,"writeVideo(v,M)");
	}
	engEvalString(ep,"close(v)");
	return 0;
}
