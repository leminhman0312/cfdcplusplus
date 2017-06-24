#include <cmath>
#include <stdio.h>
#include "engine.h"
#include <string.h>
const int nx = 10; //# in x
const int ny = 10; //# in y
const int niter = 100; //# of iterations 
#define BUFSIZE 256


int main(){
	double dx = 2/double((nx-1));
	double dy = 1/double((ny-1));
	int xmax = 2; int xmin = 0;
	int ymax = 1; int ymin = 0;
	double p[nx][ny];
	double pn[nx][ny];
	double x[nx];
	double y[ny];
	
	//populate p with zeros 

	for (int xi = 0; xi <nx;xi++){
		for (int yi = 0;yi<ny;yi++){
			p[xi][yi] = 0;			
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

	//initial condition
	
	for (int yrange = 0; yrange<ny;yrange++){
		p[nx-1][yrange] = y[yrange];
	}


	/* NO VIDEOS
	
	//WORKING LOOP!!!


	

	for (int iter = 0; iter<niter; iter++){
		//copy values
		for (int xiter= 0; xiter<nx;xiter++){
			for (int yiter = 0; yiter<ny;yiter++){
				pn[xiter][yiter]=p[xiter][yiter];
			}
		}

		//main loop
		for (int i = 1; i<nx-1;i++){
			for (int j = 1; j<ny-1;j++){
				p[i][j] = ((pow(dy,2)*(pn[i+1][j]+pn[i-1][j]))+(pow(dx,2)*(p[i][j+1]+pn[i][j-1])))/(2*(pow(dx,2)+pow(dy,2)));
			}
		}

		//boundary conditions;
		
		for (int xrange = 0; xrange<nx;xrange++){
			p[xrange][0] = p[xrange][1];
			p[xrange][ny-1] = p[xrange][ny-2];
		}

	}
	
	
	
	//Testing matrix

	for (int xs = 0; xs<nx;xs++){
		for (int ys =0; ys<ny;ys++){
			printf("%1.3f\t",p[xs][ys]);
			//printf("%f\n", x[xs]);
		}
		printf("\n");
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
	engEvalString(ep,"title('2D LAPLACE');");
	engEvalString(ep,"saveas(f,'2D_LAPLACE.png');");
	fgetc(stdin);

	*/ 

	//VIDEO PARTS

	Engine *ep = engOpen(NULL);
	mxArray *x_pointer = NULL;
	mxArray *y_pointer = NULL;
	mxArray *p_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];
	//FOR VIDEO
	engEvalString(ep,"v = VideoWriter('2D_LAPLACE.avi');");
	engEvalString(ep,"open(v)");

	//CALCULATIONS with VIDEO

	for (int iter = 0; iter<niter; iter++){
		//copy values
		for (int xiter= 0; xiter<nx;xiter++){
			for (int yiter = 0; yiter<ny;yiter++){
				pn[xiter][yiter]=p[xiter][yiter];
			}
		}

		//main loop
		for (int i = 1; i<nx-1;i++){
			for (int j = 1; j<ny-1;j++){
				p[i][j] = ((pow(dy,2)*(pn[i+1][j]+pn[i-1][j]))+(pow(dx,2)*(p[i][j+1]+pn[i][j-1])))/(2*(pow(dx,2)+pow(dy,2)));
			}
		}

		//boundary conditions;
		
		for (int xrange = 0; xrange<nx;xrange++){
			p[xrange][0] = p[xrange][1];
			p[xrange][ny-1] = p[xrange][ny-2];
		}

		//PLOTTING

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
		engEvalString(ep,"view(-45,60);");
		engEvalString(ep,"shading interp;");
		engEvalString(ep,"colormap(jet);");
		engEvalString(ep,"colorbar('southoutside');");
		engEvalString(ep,"xlabel('X VALUES');");
		engEvalString(ep,"ylabel('Y VALUES');");
		engEvalString(ep,"zlabel('PRESSURE');");
		engEvalString(ep,"title('2D LAPLACE');");
		engEvalString(ep,"M = getframe(gcf);");
		engEvalString(ep,"writeVideo(v,M)");
	}
	engEvalString(ep,"close(v)");	
	return 0;
}