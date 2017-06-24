//Solving the 2D Burgers' Equation 

// dV/dt + udV/dx + vdV/dy = nu*(d^2(u)dx^2 + d^2(v)/dx^2)

#include <iostream>
#include "engine.h"
#include <string.h>
#include <cmath>

using namespace std;

#define BUFSIZE 256


double centralDiff(double front, double middle, double end){
	return front-2*middle+end;
}

double backDiff(double first, double second){
	return first - second;
}


int main(){

	//SET UP VARIABLES
	int nx = 41;
	int ny = 41;
	int nt = 50;
	double nu = 0.01;
	double dx = 2/double((nx-1));
	double dy = 2/double((ny-1));
	double sigma = 0.009;
	double dt = (sigma*dx*dy)/(nu);
	double x[nx], y[ny]; //spatial vectors X and Y
	double u[ny][nx]; //Real solution at n+1
	double u0[ny][nx]; //Initial U, storage
	double v[ny][nx]; //Real solution at n+1
	double v0[ny][nx]; //Initial V, storage

	double constX1 = ((nu*dt)/(pow(dx,2)));
	double constY1 = ((nu*dt)/(pow(dy,2)));
	double constX2 = (dt/dx);
	double constY2 = (dt/dy);

	

	//SET UP V, V0 and U,U0 to be all 1s
	for (int i = 0; i<nx;i++){
		for (int j = 0; j<ny;j++){
			u[j][i] =1.0;
			u0[j][i]=1.0;
			v[j][i] =1.0;
			v0[j][i] =1.0;
		}
	}

	


	//SET UP INITIAL CONDITIONS 
	for (int j = 0; j<ny; j++){
		y[j] = j*dy;
		for (int i = 0; i <nx;i++){
			x[i] = i*dx;
			if((x[i]>= 0.5 && x[i]<1.0) && (y[j] >=0.5 && y[j] < 1.0)){
				u0[j][i] = 2.0;
				v0[j][i] = 2.0;
			}
			else{
				u0[j][i] = 1.0;
				v0[j][i] = 1.0;
			}
		}
	}	

	/*NO VIDEOS

	//CALCULATIONS
	for (int n = 0; n<nt;n++){
		//COPY FOR U
		for (int j = 1; j<ny;j++){
			for (int i = 1; i<nx;i++){
				u[j][i] = u0[j][i]; //copy
			}
		}

		//COPY FOR V
		for (int k = 1; k<ny;k++){
			for (int l = 1; l<nx;l++){
				v[k][l] = v0[k][l]; //copy
			}
		}


		//FINITE DIFFERENCE LOOP (not 1st and last)

		for (int j = 1; j<ny;j++){
			for (int i = 1; i<nx;i++){
				u0[j][i] = ((constX1)*(centralDiff(u[j][i+1],u[j][i],u[j][i-1])))+((constY1)*(centralDiff(u[j+1][i],u[j][i],u[j-1][i])))-(((constX2)*(u[j][i])*(backDiff(u[j][i],u[j][i-1]))))-((constY2)*(v[j][i])*(backDiff(u[j][i],u[j-1][i]))) + u[j][i];
				v0[j][i] = ((constX1)*(centralDiff(v[j][i+1],v[j][i],v[j][i-1])))+((constY1)*(centralDiff(v[j+1][i],v[j][i],v[j-1][i])))-((constX2)*(u[j][i])*(backDiff(v[j][i],v[j][i-1])))-((constY2)*(v[j][i])*(backDiff(v[j][i],v[j-1][i]))) + v[j][i];
			}
		}

	}



	
	//MATLAB PLOTTING	
	
	Engine *ep;
	mxArray *x_pointer = NULL;
	mxArray *y_pointer = NULL;
	mxArray *u_pointer = NULL;
	mxArray *u0_pointer = NULL;
	mxArray *v_pointer = NULL;
	mxArray *v0_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//Test if Matlab is opened or not
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return EXIT_FAILURE;
	}

	//PUTTING VARIABLES



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

	v0_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
	memcpy((void *)mxGetPr(v0_pointer), (void *)v0, sizeof(v0));
	engPutVariable(ep, "V0_matlab", v0_pointer);

	v_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
	memcpy((void *)mxGetPr(v_pointer), (void *)v, sizeof(v));
	engPutVariable(ep, "V_matlab", v_pointer);

	engEvalString(ep,"[X Y] = meshgrid(X_matlab,Y_matlab);");
	engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
	engEvalString(ep,"k = surf(X_matlab,Y_matlab,U_matlab);");
	engEvalString(ep,"shading interp;");
	engEvalString(ep,"colormap(jet);");
	engEvalString(ep,"zlim([1 2]);");
	engEvalString(ep,"colorbar('southoutside');");
	engEvalString(ep,"xlabel('X');");
	engEvalString(ep,"ylabel('Y');");
	engEvalString(ep,"zlabel('V = f(u,v)');");
	engEvalString(ep,"title('2D BURGERS EQUATION');");
	engEvalString(ep,"saveas(k,'2D BURGERS EQUATION.png');");
	return 0;

	*/ 


	//VIDEO PARTS

	Engine *ep = engOpen(NULL);
	mxArray *x_pointer = NULL;
	mxArray *y_pointer = NULL;
	mxArray *u_pointer = NULL;
	mxArray *u0_pointer = NULL;
	mxArray *v_pointer = NULL;
	mxArray *v0_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//OPEN VIDEOS
	engEvalString(ep,"v = VideoWriter('BURGERS.avi');");
	engEvalString(ep,"open(v)");

	//engEvalString(ep,"zlim([1 2]);");
		
		


	for (int n = 0; n<nt;n++){
		//COPY FOR U
		for (int j = 1; j<ny;j++){
			for (int i = 1; i<nx;i++){
				u[j][i] = u0[j][i]; //copy
			}
		}

		//COPY FOR V
		for (int k = 1; k<ny;k++){
			for (int l = 1; l<nx;l++){
				v[k][l] = v0[k][l]; //copy
			}
		}


		//FINITE DIFFERENCE LOOP (not 1st and last)

		for (int j = 1; j<ny;j++){
			for (int i = 1; i<nx;i++){
				u0[j][i] = ((constX1)*(centralDiff(u[j][i+1],u[j][i],u[j][i-1])))+((constY1)*(centralDiff(u[j+1][i],u[j][i],u[j-1][i])))-(((constX2)*(u[j][i])*(backDiff(u[j][i],u[j][i-1]))))-((constY2)*(v[j][i])*(backDiff(u[j][i],u[j-1][i]))) + u[j][i];
				v0[j][i] = ((constX1)*(centralDiff(v[j][i+1],v[j][i],v[j][i-1])))+((constY1)*(centralDiff(v[j+1][i],v[j][i],v[j-1][i])))-((constX2)*(u[j][i])*(backDiff(v[j][i],v[j][i-1])))-((constY2)*(v[j][i])*(backDiff(v[j][i],v[j-1][i]))) + v[j][i];
			}
		}

		//PLOTTING with VIDEOS

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

		engEvalString(ep,"[X Y] = meshgrid(X_matlab,Y_matlab);");
		engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
		engEvalString(ep,"k = surf(X_matlab,Y_matlab,U_matlab);");
		engEvalString(ep,"shading interp;");
		engEvalString(ep,"colormap(jet);");
		engEvalString(ep,"zlim([1 2]);");
		engEvalString(ep,"colorbar('southoutside');");
		engEvalString(ep,"xlabel('X');");
		engEvalString(ep,"ylabel('Y');");
		engEvalString(ep,"zlabel('V = f(u,v)');");
		engEvalString(ep,"title('2D BURGERS EQUATION');");		
		engEvalString(ep,"M = getframe(gcf);");
		engEvalString(ep,"writeVideo(v,M)");
	}
	engEvalString(ep,"close(v)");
	return 0;

}
	
	

