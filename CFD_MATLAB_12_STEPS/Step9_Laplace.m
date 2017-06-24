% Step 9 Laplace Equation
% Manuel Ramsaier

%Step9_Laplace
clc
clear all
nx=5; ny=5; nit=100;
dx = 2/(nx-1); dy=1/(ny-1);
x=0:dx:2; y=0:dy:1;
p=zeros(nx,ny);
p(nx,:)=y;
[~,~]=meshgrid(y,x);




for iit = 1:nit
    pd=p;
    for i=2:nx-1	
    for j=2:ny-1
    p(i,j) = ((pd(i+1,j)+pd(i-1,j))*dy^2+(pd(i,j-1)+pd(i,j+1))*dx^2 )/(dx^2+dy^2)/2;
    end
    end
    p(2:nx-1,1)=p(2:nx-1,2);
    p(2:nx-1,ny)=p(2:nx-1,ny-1);

end

%surf(x,y,p)

