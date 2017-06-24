//SOLVING POISSON EQUATION 

// del^2(P) = b

#include <cmath>
#include <stdio.h>
#include "engine.h"
#include <string.h>
const int nx = 101; //# in x
const int ny = 101; //# in y
const int nt = 100; //# of iterations 
#define BUFSIZE 256


int main(){
	int xmax = 2; int xmin = 0;
	int ymax = 1; int ymin = 0;
	double dx = (xmax-xmin)/double((nx-1));
	double dy = (ymax-ymin)/double((ny-1));
	double p[nx][ny];
	double pn[nx][ny];
	double x[nx];
	double y[ny];
	double b[nx][ny];

	int threeforthnx = (0.75 * (nx-1));
	int threeforthny = (0.75 * (ny-1));
	int forthnx = (0.25 * (nx-1));
	int forthny = (0.25 * (ny-1));

	//populate p,pn, with zeros 

	for (int xi = 0; xi <nx;xi++){
		for (int yi = 0;yi<ny;yi++){
			p[xi][yi] = 0;	
			pn[xi][yi] = 0;		
			b[xi][yi] = 0;
		}
	}

	//populate x and y

	//X
	for (int xnum = 0; xnum <nx; xnum++){
		x[0] = 0;
		x[xnum+1] = x[xnum] + dx;
	}

	//Y
	for (int ynum = 0; ynum <ny; ynum++){
		y[0] = 0;
		y[ynum+1] = y[ynum] + dy;
	}


	//POPULATE INITIAL CONDITION b

	b[threeforthnx][threeforthny] = -100;
	b[forthnx][forthny] = 100;

	/* NO VIDEOS

	//Calculations
	
	for (int time = 0; time <nt; time++){
		//copy loops
		for (int xiter= 0; xiter<nx;xiter++){
			for (int yiter = 0; yiter<ny;yiter++){
				pn[xiter][yiter]=p[xiter][yiter];
			}
		}

		//main loop
		for (int i = 1; i<nx-1;i++){
			for (int j = 1; j<ny-1;j++){
				p[i][j] = ((pow(dy,2)*(pn[i+1][j]+pn[i-1][j]))+(pow(dx,2)*(p[i][j+1]+pn[i][j-1]))-(b[i][j]*pow(dx,2)*pow(dy,2)))/(2*(pow(dx,2)+pow(dy,2)));
			}
		}

		for (int xrange = 0; xrange<nx;xrange++){
			p[xrange][0] = p[xrange][1];
			p[xrange][ny-1] = p[xrange][ny-2];
		}
	}

	//MATLAB PLOTTING
	
	
	
	Engine *ep;
	mxArray *x_pointer = NULL;
	mxArray *y_pointer = NULL;
	mxArray *p_pointer = NULL, *result = NULL;
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

	p_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
	memcpy((void *)mxGetPr(p_pointer), (void *)p, sizeof(p));
	engPutVariable(ep, "P_matlab", p_pointer);
	

	engEvalString(ep,"[~,~] = meshgrid(X_matlab,X_matlab);");
	//engEvalString(ep,"set(gcf, 'Visible', 'off');");
	engEvalString(ep,"f = surf(X_matlab,Y_matlab,P_matlab);");
	engEvalString(ep,"axis([0 2 0 1]);");
	engEvalString(ep,"view(-45,60);");
	engEvalString(ep,"shading interp;");
	engEvalString(ep,"colormap(jet);");
	engEvalString(ep,"colorbar('southoutside');");
	engEvalString(ep,"xlabel('X VALUES');");
	engEvalString(ep,"ylabel('Y VALUES');");
	engEvalString(ep,"zlabel('PRESSURE');");
	engEvalString(ep,"title('2D POISSON');");
	engEvalString(ep,"saveas(f,'2D_POISSON.png');");
	fgetc(stdin);
	
	*/ 

	//WITH VIDEOS

	Engine *ep = engOpen(NULL);
	mxArray *x_pointer = NULL;
	mxArray *y_pointer = NULL;
	mxArray *p_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//FOR VIDEO

	engEvalString(ep,"v = VideoWriter('2D_POISSON.avi');");
	engEvalString(ep,"open(v)");	

	//Calculations with VIDEO

	for (int time = 0; time <nt; time++){
		//copy loops
		for (int xiter= 0; xiter<nx;xiter++){
			for (int yiter = 0; yiter<ny;yiter++){
				pn[xiter][yiter]=p[xiter][yiter];
			}
		}

		//main loop
		for (int i = 1; i<nx-1;i++){
			for (int j = 1; j<ny-1;j++){
				p[i][j] = ((pow(dy,2)*(pn[i+1][j]+pn[i-1][j]))+(pow(dx,2)*(p[i][j+1]+pn[i][j-1]))-(b[i][j]*pow(dx,2)*pow(dy,2)))/(2*(pow(dx,2)+pow(dy,2)));
			}
		}

		for (int xrange = 0; xrange<nx;xrange++){
			p[xrange][0] = p[xrange][1];
			p[xrange][ny-1] = p[xrange][ny-2];
		}

		x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
		memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
		engPutVariable(ep, "X_matlab", x_pointer);
		y_pointer = mxCreateDoubleMatrix(1,ny,mxREAL);
		memcpy((void *)mxGetPr(y_pointer), (void *)y, sizeof(y));
		engPutVariable(ep, "Y_matlab", y_pointer);
		p_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
		memcpy((void *)mxGetPr(p_pointer), (void *)p, sizeof(p));
		engPutVariable(ep, "P_matlab", p_pointer);
		engEvalString(ep,"[~,~] = meshgrid(X_matlab,X_matlab);");
		engEvalString(ep,"set(gcf, 'Visible', 'off');");
		engEvalString(ep,"f = surf(X_matlab,Y_matlab,P_matlab);");
		engEvalString(ep,"axis([0 2 0 1]);");
		engEvalString(ep,"view(90,0);");
		engEvalString(ep,"shading interp;");
		engEvalString(ep,"colormap(jet);");
		engEvalString(ep,"colorbar('southoutside');");
		engEvalString(ep,"xlabel('X VALUES');");
		engEvalString(ep,"ylabel('Y VALUES');");
		engEvalString(ep,"zlabel('PRESSURE');");
		engEvalString(ep,"title('2D POISSON');");
		engEvalString(ep,"M = getframe(gcf);");
		engEvalString(ep,"writeVideo(v,M)");
	}
	engEvalString(ep,"close(v)");
	return 0;
}