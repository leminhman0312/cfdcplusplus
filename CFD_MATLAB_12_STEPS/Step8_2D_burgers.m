%%%variable declarations
clear all
nx = 41;
ny = 41;
nt = 120;
c=1;
dx = 2.0/(nx-1);
dy = 2.0/(ny-1);
sigma = .009;
nu=.01;
dt = sigma*dx*dy/nu;

x = linspace(0,2,nx);
y = linspace(0,2,ny);

u = ones(ny,nx); %%create a 1xn vector of 1's
v = ones(ny,nx);
un=ones(ny,nx);
vn=ones(ny,nx);
comb=ones(ny,nx);

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
        u(i,j)=u(i,j)- (dt/dx) * u(i,j)*(u(i,j) -u(i-1,j)) - (dt/dy) * v(i,j)*(u(i,j)-u(i,j-1)) + (nu*dt/dx^2) *(u(i-1,j)-2*u(i,j)+u(i-1,j)) + (nu*dt/dy^2) * (u(i,j+1)-2*u(i,j)+u(i,j-1));
        v(i,j)=v(i,j)- (dt/dx) * u(i,j)*(v(i,j) -v(i-1,j)) - (dt/dy) * v(i,j)*(v(i,j)-v(i,j-1)) + (nu*dt/dx^2) *( v(i-1,j)-2*v(i,j)+v(i-1,j)) + (nu*dt/dy^2) * (v(i,j+1)-2*v(i,j)+v(i,j-1));
            u(1:ny,1)=1;
            u(1,1:nx)=1;
            u(ny,1:nx)=1;
            u(1:nx,ny)=1;
            v(1:ny,1)=1;
            v(1,1:nx)=1;
            v(ny,1:nx)=1;
            v(1:nx,ny)=1;
        end
    end
  
   %pause(0.1)
end

for n = 1:nt+1
     surf(x,y,u)
     drawnow()
end



