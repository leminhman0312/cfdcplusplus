%%%variable declarations
clear all
nx = 101;
ny = 101;
nt = 80;
c = 1;
dx = 2.0/(nx-1);
dy = 2.0/(ny-1);
sigma = .2;
dt = sigma*dx;

x = linspace(0,2,nx);
y = linspace(0,2,ny);

u = ones(ny,nx); %%create a 1xn vector of 1's
v=ones(ny,nx);
un=ones(ny,nx);
vn=ones(ny,nx);

un = ones(ny,nx); %%

%%%Assign initial conditions

u(.5/dy:1/dy+1,.5/dx:1/dx+1)=2; %%set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
v(.5/dy:1/dy+1,.5/dx:1/dx+1)=2; %%set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

%%%Plot Initial Condition
 %%the figsize parameter can be used to produce different sized images
                  
[X, Y] = meshgrid(x,y);                            


for n=1:nt+1
    un=u;
    vn=v;
    
    for i=2:(ny-1)
        for j=2:(nx-1)
        
            u(i,j)=un(i,j)-(un(i,j)*dt/dx*(un(i,j)-un(i-1,j)))-vn(i,j)*dt/dy*(un(i,j)-un(i,j-1));
            v(i,j)=vn(i,j)-(un(i,j)*dt/dx*(vn(i,j)-vn(i-1,j)))-vn(i,j)*dt/dy*(vn(i,j)-vn(i,j-1));
            u(1:ny,1)=1;
            u(1,1:nx)=1;
            u(ny,1:nx)=1;
            u(1:ny,nx)=1;
            
            v(1:ny,1)=1;
            v(1,1:nx)=1;
            v(ny,1:nx)=1;
            v(1:ny,nx)=1;
            
            
        end
    end
   surf(x,y,u)
   pause(0.1)
end

