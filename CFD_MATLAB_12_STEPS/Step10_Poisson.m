%Step 10 :  Poisson Equation

% 

clear all
nx=20; ny=20; nt=100;
xmin=0; xmax=2; 
ymin=0; ymax=1;
dx = (xmax-xmin)/(nx-1); 
dy=(ymax-ymin)/(ny-1);


% Init
p=zeros(nx,ny);
pd=zeros(nx,ny);
b=zeros(nx,ny);
x=xmin:dx:xmax; 
y=ymin:dy:ymax;

% Source

b(floor(nx/4),floor(ny/4))=100;
b(floor(3*nx/4),floor(3*ny/4))=-100;

for nt = 1:nt
    pd=p;
    for i=2:nx-1
    for j=2:ny-1
    p(i,j) = ((pd(i+1,j)+pd(i-1,j))*dy^2+ (pd(i,j-1)+pd(i,j+1))*dx^2 -b(i,j)*dx^2*dy^2 )/(dx^2+dy^2)/2;
    end
    end
    p(2:nx-1,1)=p(2:nx-1,2);
    p(2:nx-1,ny)=p(2:nx-1,ny-1);

end

surf(x,y,p)
