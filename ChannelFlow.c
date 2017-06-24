//SOLVING Channel Flow in 2D
// U is in Y dir, V is in X dir

#include <stdio.h>
#include <math.h>
#include "engine.h"
#include <string.h>
#define BUFSIZE 256
double x_vector(double x[], int nx, double dx);
double secondOrderCentral(double front, double middle, double last, double deltaSecond);
double forwardbackDiff(double first, double second, double deltaF);
double centralDiffX(double firstX, double secondX, double deltaX);
double centralDiffY(double firstY, double secondY, double deltaY);
double y_vector(double y[], int ny, double dy);
double SquareBracketPoisson(int nx, int ny, double b[][ny], double rho, double dt, double u[][ny], double v[][ny], double dx, double dy,double onebydt);
double PressurePoisson(int nx, int ny,  double pn[][ny],double p[][ny], double b[][ny], double dx, double dy, int nit,double dxsquare, double dysquare);
double CopyLoop(int nx, int ny, double real[][ny], double copy[][ny]);
double ChannelFlow(int nt, int nx, int ny, double dx, double dy, double dt, double rho, double nu, double b[][ny], double pn[][ny],int nit,double vn[][ny], double un[][ny],double onebydt,double dxsquare, double dysquare,double Fsource);


int main(){

	//DECLARE VARIABLES
	int nx = 41; //step in x
	int ny = 41; //step in y
	int nt = 50; //step in time
	int nit = 50; //iterations
	int c = 1.; 
	double dx = 2/(double)(nx-1);
	double dy = 2/(double)(ny-1);
	double x[nx]; //x vector
	double y[ny]; //y vector
	double rho = 1.; //density
	double nu = 0.1; //viscosity

	double dt = 0.001; //time 
	double onebydt = 1/dt;
	double dxsquare = pow(dx,2);
	double dysquare = pow(dy,2);
	double Fsource = 1.0;
	double u[nx][ny]; //U component
	double v[nx][ny]; //V component
	double un[nx][ny]; //holder for U
	double vn[nx][ny]; //holder for V
	double b[nx][ny]; //square brackets of Pressure Poisson Equation
	double p[nx][ny]; //Pressure field
	double pn[nx][ny]; //Pressure field


	//Setting zero matrices (all rows, columns)

	for (int x = 0; x<nx;x++){
		for (int y = 0; y<ny;y++){
			u[x][y] = 0;
			v[x][y] = 0;
			b[x][y] = 0;
			p[x][y] = 0;
			pn[x][y] = 0;
			un[x][y] = 0;
			vn[x][y] = 0;
		}
	}
	


	//CALCULATE CAVITY FLOW in 2D
	
	ChannelFlow(nt,nx,ny,dx,dy,dt,rho,nu,b,pn,nit,vn,un,onebydt,dxsquare,dysquare,Fsource);
	

	return 0;

	
}


//FUNCTIONS

double centralDiffX(double firstX, double secondX, double deltaX){
	return (firstX - secondX)/(2*deltaX);
}

double centralDiffY(double firstY, double secondY, double deltaY){
	return (firstY - secondY)/(2*deltaY);
}

double forwardbackDiff(double first, double second, double deltaF){
	return (first - second)/(deltaF);
}

double secondOrderCentral(double front, double middle, double last, double deltaSecond){
	return (front -(2*middle) + last)/(pow(deltaSecond,2));
}

double x_vector(double x[], int nx, double dx){
	for (int i = 0; i<nx; i++){
		x[0] = 0;
		x[i+1] = x[i] + dx;
	}
}

double y_vector(double y[], int ny, double dy){
	for (int j = 0; j<ny; j++){
		y[0] = 0;
		y[j+1] = y[j] + dy;
	}
}

double SquareBracketPoisson(int nx, int ny, double b[][ny], double rho, double dt, double u[][ny], double v[][ny], double dx, double dy, double onebydt){
	//looping, excluding first and last elements

	for (int i = 1; i<nx-1;i++){
		for (int j = 1; j<ny-1;j++){
			b[i][j] = rho*((onebydt*(centralDiffX(u[i+1][j],u[i-1][j],dx)+centralDiffY(v[i][j+1],v[i][j-1],dy)))-(pow(centralDiffX(u[i+1][j],u[i-1][j],dx),2))+(2*centralDiffY(u[i][j+1],u[i][j-1],dy)*centralDiffX(v[i+1][j],v[i-1][j],dx))+(pow(centralDiffY(v[i][j+1],v[i][j-1],dy),2)));
			}
		}

		
		//Periodic Conditions

		for (int i = 1; i<nx-1;i++){
			for (int j = 1; j<ny-1;j++){
				//X = 0
				b[0][j] = rho*((onebydt*(centralDiffX(u[1][j],u[nx-1][j],dx)+centralDiffY(v[i][j+1],v[i][j-1],dy)))-(pow(centralDiffX(u[1][j],u[nx-1][j],dx),2))+(2*centralDiffY(u[i][j+1],u[i][j-1],dy)*centralDiffX(v[i+1][j],v[i-1][j],dx))+(pow(centralDiffY(v[i][j+1],v[i][j-1],dy),2)));
				//X = 2
				b[nx-1][j] = rho*((onebydt*(centralDiffX(u[0][j],u[nx-2][j],dx)+centralDiffY(v[i][j+1],v[i][j-1],dy)))-(pow(centralDiffX(u[0][j],u[nx-2][j],dx),2))+(2*centralDiffY(u[i][j+1],u[i][j-1],dy)*centralDiffX(v[i+1][j],v[i-1][j],dx))+(pow(centralDiffY(v[i][j+1],v[i][j-1],dy),2)));
			}
		}	
}





double CopyLoop(int nx, int ny, double real[][ny], double copy[][ny]){
	for (int xiter = 0; xiter<nx;xiter++){
			for (int yiter = 0; yiter<ny; yiter++){
				copy[xiter][yiter] = real[xiter][yiter];
			}
		}
}


double PressurePoisson(int nx, int ny, double pn[][ny], double p[][ny],double b[][ny], double dx, double dy, int nit, double dxsquare, double dysquare){
	for (int iit = 0; iit <=nit; iit++){

		//COPY P to PN AS PLACEHOLDER
		CopyLoop(nx,ny,p,pn);
		
		//CALCULATE P
		for (int i = 1; i<nx-1;i++){
			for (int j = 1; j<ny-1;j++){
				p[i][j] = ((((pn[i+1][j]+pn[i-1][j])*(dysquare))+((pn[i][j+1]+pn[i][j-1])*(dxsquare)))/(2*(dxsquare+dysquare)))-((dxsquare*dysquare)/(2*(dxsquare+dysquare)))*b[i][j];
			}
		}

		//PERIODIC CONDITIONS FOR P

		for (int i = 1; i<nx-1;i++){
			for (int j = 1; j<ny-1;j++){
				//for X = 0 for P
				p[0][j] = ((((pn[1][j]+pn[nx-1][j])*(dysquare))+((pn[0][j+1]+pn[0][j-1])*(dxsquare)))/(2*(dxsquare+dysquare)))-((dxsquare*dysquare)/(2*(dxsquare+dysquare)))*b[i][j];
			}
		}

		for (int i =1; i<nx-1;i++){
			for (int j = 1; j<ny-1; j++){
				//for X = 2 for P
				p[nx-1][j] = ((((pn[0][j]+pn[nx-2][j])*(dysquare))+((pn[nx-1][j+1]+pn[nx-1][j-1])*(dxsquare)))/(2*(dxsquare+dysquare)))-((dxsquare*dysquare)/(2*(dxsquare+dysquare)))*b[i][j];

			}
		}


		//BOUNDARY CONDITIONS FOR P
		for (int rangep = 0; rangep<nx; rangep++){
			p[rangep][nx-1] = p[rangep][nx-2];
			//p[0][rangep] = p[1][rangep];
			p[rangep][0] = p[rangep][1];
			//p[ny-1][rangep] = 0.0;
		}
	}
}

double ChannelFlow(int nt, int nx, int ny, double dx, double dy, double dt, double rho, double nu, double b[][ny], double pn[][ny],int nit,double vn[][ny], double un[][ny],double onebydt,double dxsquare, double dysquare, double Fsource){
	
	double u[nx][ny];
	double v[nx][ny];
	double p[nx][ny];
	double x[nx];
	double y[ny];
	double col_u[nx];
	double col_un[nx];
	double sum_col_u = 0;
	double sum_col_un = 0;
	double udiff = 1;
	double stepcount = 0;
	


	//POPULATE X AND Y VECTORS
	for (int i = 0; i<nx; i++){
		x[0] = 0;
		x[i+1] = x[i] + dx;
		y[0] = 0;
		y[i+1] = y[i] + dy;
	}

	//MATLAB DECLARE POINTERS
	Engine *ep = engOpen(NULL);
	mxArray *p_pointer = NULL;
	mxArray *x_pointer = NULL;
	mxArray *u_pointer = NULL;
	mxArray *v_pointer = NULL;
	mxArray *y_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//Open Videos
	engEvalString(ep,"v = VideoWriter('CHANNEL_2D.avi');");
	//engEvalString(ep,"v.FrameRate = 30");
	engEvalString(ep,"open(v)");


	//INSTEAD OF DOING THIS FOR LOOP
	//for (int time = 0; time<nt;time++){
	// We can do a while loop instead

	while(udiff > 0.001){

		//CALL FUNCTION TO CALCULATE B
		SquareBracketPoisson(nx,ny,b,rho,dt,u,v,dx,dy,onebydt);
		//CALL FUNCTION TO CALCULATE P
		PressurePoisson(nx, ny, pn, p, b, dx, dy,nit,dxsquare,dysquare);
		//COPY U, V to UN, VN AS PLACEHOLDERS
		CopyLoop(nx,ny,u,un);
		CopyLoop(nx,ny,v,vn);
		//CALCULATE VELOCITY FIELD
		for (int i =1; i<nx-1;i++){
			for (int j = 1; j<ny-1; j++){
				u[i][j] = un[i][j]-(un[i][j]*dt*(forwardbackDiff(un[i][j],un[i-1][j],dx)))-(vn[i][j]*dt*(forwardbackDiff(un[i][j],un[i][j-1],dy)))-((dt/(rho))*(centralDiffX(p[i+1][j],p[i-1][j],dx)))+((nu*dt)*(secondOrderCentral(un[i+1][j],un[i][j],un[i-1][j],dx)))+((nu*dt)*(secondOrderCentral(un[i][j+1],un[i][j],un[i][j-1],dy)))+(Fsource*dt);

				v[i][j] = vn[i][j]-(un[i][j]*dt*(forwardbackDiff(vn[i][j],vn[i-1][j],dx)))-(vn[i][j]*dt*(forwardbackDiff(vn[i][j],vn[i][j-1],dy)))-((dt/(rho))*(centralDiffX(p[i][j+1],p[i][j-1],dy)))+((nu*dt)*(secondOrderCentral(vn[i+1][j],vn[i][j],vn[i-1][j],dx)))+((nu*dt)*(secondOrderCentral(vn[i][j+1],vn[i][j],vn[i][j-1],dy)));
			}
		}



		//PERIODIC CONDITIONS FOR U and V

		for (int i =1; i<nx-1;i++){
			for (int j = 1; j<ny-1; j++){
				//X = 0 for U

				u[0][j] = un[0][j]-(un[0][j]*dt*(forwardbackDiff(un[0][j],un[nx-1][j],dx)))-(vn[0][j]*dt*(forwardbackDiff(un[0][j],un[0][j-1],dy)))-((dt/(rho))*(centralDiffX(p[1][j],p[nx-1][j],dx)))+((nu*dt)*(secondOrderCentral(un[1][j],un[0][j],un[nx-1][j],dx)))+((nu*dt)*(secondOrderCentral(un[0][j+1],un[0][j],un[0][j-1],dy)))+(Fsource*dt);
			}
		}

		for (int i =1; i<nx-1;i++){
			for (int j = 1; j<ny-1; j++){
				//X = 0 for V

				v[0][j] = vn[0][j]-(un[0][j]*dt*(forwardbackDiff(vn[0][j],vn[nx-1][j],dx)))-(vn[0][j]*dt*(forwardbackDiff(vn[0][j],vn[0][j-1],dy)))-((dt/(rho))*(centralDiffX(p[0][j+1],p[0][j-1],dy)))+((nu*dt)*(secondOrderCentral(vn[1][j],vn[0][j],vn[nx-1][j],dx)))+((nu*dt)*(secondOrderCentral(vn[0][j+1],vn[0][j],vn[0][j-1],dy)));
			}
		}

		for (int i =1; i<nx-1;i++){
			for (int j = 1; j<ny-1; j++){

				//X = 2 for U;

				u[nx-1][j] = un[nx-1][j]-(un[nx-1][j]*dt*(forwardbackDiff(un[nx-1][j],un[nx-2][j],dx)))-(vn[nx-1][j]*dt*(forwardbackDiff(un[nx-1][j],un[nx-1][j-1],dy)))-((dt/(rho))*(centralDiffX(p[0][j],p[nx-2][j],dx)))+((nu*dt)*(secondOrderCentral(un[0][j],un[nx-1][j],un[nx-2][j],dx)))+((nu*dt)*(secondOrderCentral(un[nx-1][j+1],un[nx-1][j],un[nx-1][j-1],dy)))+(Fsource*dt);
			}
		}

		for (int i =1; i<nx-1;i++){
			for (int j = 1; j<ny-1; j++){

				//X = 2 for V

				v[nx-1][j] = vn[nx-1][j]-(un[nx-1][j]*dt*(forwardbackDiff(vn[nx-1][j],vn[nx-2][j],dx)))-(vn[nx-1][j]*dt*(forwardbackDiff(vn[nx-1][j],vn[nx-1][j-1],dy)))-((dt/(rho))*(centralDiffX(p[nx-1][j+1],p[nx-1][j-1],dy)))+((nu*dt)*(secondOrderCentral(vn[0][j],vn[nx-1][j],vn[nx-2][j],dx)))+((nu*dt)*(secondOrderCentral(vn[nx-1][j+1],vn[nx-1][j],vn[nx-1][j-1],dy)));
			}
		}


		//SETTING UP VELOCITY BOUNDARY CONDITIONS
		//(U goes with Y (j), V goes with X(i))

		//U, V = 0 at y = 0, y = 2
		for (int xbound = 0; xbound <nx; xbound++){
			for (int ybound = 0; ybound <ny; ybound++){
				u[xbound][0] = 0;
				v[xbound][0] = 0;
				u[xbound][ny-1] = 0;
				v[xbound][ny-1] = 0;
			}
		}

		//Calculate SUM(U) and SUM(UN)

		//Sum of columns first

		for (int j = 0; j<ny;j++){
			col_u[j] = 0;
			col_un[j] = 0;
			for (int i = 0; i<nx;i++){
				col_u[j] = col_u[j] + u[i][j];
				col_un[j] = col_un[j] + un[i][j];
			}
		}

		//Sum of sum of columns

		for (int step = 0; step <nx;step++){
			sum_col_u = sum_col_u+col_u[step];
			sum_col_un = sum_col_un+col_un[step];
		}

		//Calculate udiff

		udiff = (sum_col_u-sum_col_un)/(sum_col_u);
		stepcount = stepcount + 1;

		//PUT INTO MATLAB
		x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
		memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
		engPutVariable(ep, "X_matlab", x_pointer);
		y_pointer = mxCreateDoubleMatrix(1,ny,mxREAL);
		memcpy((void *)mxGetPr(y_pointer), (void *)y, sizeof(y));
		engPutVariable(ep, "Y_matlab", y_pointer);
		p_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
		memcpy((void *)mxGetPr(p_pointer), (void *)p, sizeof(p));
		engPutVariable(ep, "P_matlab", p_pointer);
		u_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
		memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
		engPutVariable(ep, "U_matlab", u_pointer);
		v_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
		memcpy((void *)mxGetPr(v_pointer), (void *)v, sizeof(v));
		engPutVariable(ep, "V_matlab", v_pointer);

		//PLOTTING

		//MESH GRID
		engEvalString(ep,"[X Y] = meshgrid(X_matlab,Y_matlab);");
		engEvalString(ep,"set(gcf, 'Visible', 'off');"); //hidden plot figure
		engEvalString(ep,"quiver(X,Y,U_matlab.',V_matlab.',1);");
		engEvalString(ep,"pause(0.01)");
		engEvalString(ep,"xlabel('X');");
		engEvalString(ep,"ylabel('Y');");
		engEvalString(ep,"title('Chanel Flow 2D');");
		//engEvalString(ep,"axis([0 2 0 2]);");
		//engEvalString(ep,"saveas(c,'channelflow2d.png');");
		engEvalString(ep,"M = getframe(gcf);");
		engEvalString(ep,"writeVideo(v,M)");
	}
	engEvalString(ep,"close(v)");
	mxDestroyArray(x_pointer);
	mxDestroyArray(y_pointer);
	mxDestroyArray(p_pointer);
	mxDestroyArray(u_pointer);
	mxDestroyArray(v_pointer);
	engEvalString(ep,"close");
}




		










	