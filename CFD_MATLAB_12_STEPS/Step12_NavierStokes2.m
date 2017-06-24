%Step 10 :  Poisson Equation

% 

clear all
nx=21; ny=21; nt=50; nit=50;
xmin=0; xmax=2; 
ymin=0; ymax=2;
dx = (xmax-xmin)/(nx-1); 
dy=(ymax-ymin)/(ny-1);
x=linspace(0,2,nx);
y=linspace(0,2,ny);
[Y,X]=meshgrid(y,x);

rho=1;
nu=0.1;
F=1;
dt=0.01;


% Init
u=zeros(ny,nx);
v=zeros(ny,nx);
p=zeros(ny,nx);
b=zeros(ny,nx);
%Pressure Field
%Square Brackets of Poissons Equation
  

udiff=1;
stepcount=0;
while udiff >0.001
        for i=2:(nx-1)
            for j=2:(ny-1)
                b(i,j)=rho*(1/dt*((u(i+1,j)-u(i-1,j))/(2*dx)+(v(i,j+1)-v(i,j-1))/(2*dy))-((u(i+1,j)-u(i-1,j))/(2*dx))^2-2*((u(i,j+1)-u(i,j-1))/(2*dy)*(v(i+1,j)-v(i-1,j))/(2*dx))-((v(i,j+1)-v(i,j-1))/(2*dy))^2);
            end
        end
        %Periodischer Term für x=0 ( Was bei x= 0 passiert, passiert bei
        %x=2)
            for j=2:(ny-1)
                b(1,j)=rho*(1/dt*((u(2,j)-u(nx,j))/(2*dx)+(v(i,j+1)-v(i,j-1))/(2*dy))-((u(2,j)-u(nx,j))/(2*dx))^2-2*((u(i,j+1)-u(i,j-1))/(2*dy)*(v(2,j)-v(nx,j))/(2*dx))-((v(i,j+1)-v(i,j-1))/(2*dy))^2);
            end
            
        % Periodischer Term für x = 2
            for j=2:ny-1
                b(nx,j)=rho*(1/dt*((u(1,j)-u(nx-1,j))/(1*dx)+(v(i,j+1)-v(i,j-1))/(1*dy))-((u(1,j)-u(nx-1,j))/(1*dx))^1-1*((u(i,j+1)-u(i,j-1))/(1*dy)*(v(1,j)-v(nx-1,j))/(1*dx))-((v(i,j+1)-v(i,j-1))/(1*dy))^1);
            end
        
        
    for iit=1:nit+1
        pn=p;
        for i=2:(nx-1)
            for j=2:(ny-1) 
            p(i,j)=((pn(i+1,j)+pn(i-1,j))*dy^2+(pn(i,j+1)+pn(i,j-1))*dx^2)/(2*(dx^2+dy^2))-dx^2*dy^2/(2*(dx^2+dy^2))*b(i,j);
            end
        end
        %Periodischer Term für x=0
        
            for j=2:(ny-1) 
            p(1,j)=((pn(2,j)+pn(nx,j))*dy^2+(pn(1,j+1)+pn(1,j-1))*dx^2)/(2*(dx^2+dy^2))-dx^2*dy^2/(2*(dx^2+dy^2))*b(i,j);
            end
        
        %Periodischer Term für x=2
            for j=2:(ny-1) 
            p(nx,j)=((pn(1,j)+pn(nx-1,j))*dy^2+(pn(nx,j+1)+pn(nx,j-1))*dx^2)/(2*(dx^2+dy^2))-dx^2*dy^2/(2*(dx^2+dy^2))*b(i,j);
            end
        
        p(:,ny) =p(:,ny-1);	%%dp/dy = 0 at y = 2
        p(:,1) = p(:,2);		%%dp/dy = 0 at y = 0

    end
    
    un = u;
    vn = v;
    
    for j=2:nx-1
        for i=2:ny-1        
        %Velocity Field
        u(i,j) = un(i,j)-un(i,j)*dt/dx*(un(i,j)-un(i-1,j))-vn(i,j)*dt/dy*(un(i,j)-un(i,j-1))-dt/(2*rho*dx)*(p(i+1,j)-p(i-1,j))+nu*(dt/dx^2*(un(i+1,j)-2*un(i,j)+un(i-1,j))+ (dt/dy^2*(un(i,j+1)-2*un(i,j)+un(i,j-1))))+F*dt;
        v(i,j) = vn(i,j)-un(i,j)*dt/dx*(vn(i,j)-vn(i-1,j))-vn(i,j)*dt/dy*(vn(i,j)-vn(i,j-1))-dt/(2*rho*dy)*(p(i,j+1)-p(i,j-1))+nu*(dt/dx^2*(vn(i+1,j)-2*vn(i,j)+vn(i-1,j))+ (dt/dy^2*(vn(i,j+1)-2*vn(i,j)+vn(i,j-1))));
        end
    end
    
    %Periodischer Term für x =0
    for j=2:ny-1
    u(1,j) = un(1,j)-un(1,j)*dt/dx*(un(1,j)-un(nx,j))-vn(1,j)*dt/dy*(un(1,j)-un(1,j-1))-dt/(2*rho*dx)*(p(2,j)-p(nx,j))+nu*(dt/dx^2*(un(2,j)-2*un(1,j)+un(nx,j))+ (dt/dy^2*(un(1,j+1)-2*un(1,j)+un(1,j-1))))+F*dt;
        v(1,j) = vn(1,j)-un(1,j)*dt/dx*(vn(1,j)-vn(nx,j))-vn(1,j)*dt/dy*(vn(1,j)-vn(1,j-1))-dt/(2*rho*dy)*(p(1,j+1)-p(1,j-1))+nu*(dt/dx^2*(vn(2,j)-2*vn(1,j)+vn(nx,j))+ (dt/dy^2*(vn(1,j+1)-2*vn(1,j)+vn(1,j-1))));
    end
    %Periodischer Term für x=2
    for j=2:ny-1
    u(nx,j) = un(nx,j)-un(nx,j)*dt/dx*(un(nx,j)-un(nx-1,j))-vn(nx,j)*dt/dy*(un(nx,j)-un(nx,j-1))-dt/(2*rho*dx)*(p(1,j)-p(nx-1,j))+nu*(dt/dx^2*(un(1,j)-2*un(nx,j)+un(nx-1,j))+ (dt/dy^2*(un(nx,j+1)-2*un(nx,j)+un(nx,j-1))))+F*dt;
    v(nx,j) = vn(nx,j)-un(nx,j)*dt/dx*(vn(nx,j)-vn(nx-1,j))-vn(nx,j)*dt/dy*(vn(nx,j)-vn(nx,j-1))-dt/(2*rho*dy)*(p(nx,j+1)-p(nx,j-1))+nu*(dt/dx^2*(vn(1,j)-2*vn(nx,j)+vn(nx-1,j))+ (dt/dy^2*(vn(nx,j+1)-2*vn(nx,j)+vn(nx,j-1))));
    end
%Wall BC: u,v = 0 @ y = 0,2
u(:,ny)=0;
u(:,1)=0;
v(:,1)=0;
v(:,ny)=0;
udiff = (sum(sum(u))-sum(sum(un)))/sum(sum(u))
stepcount=stepcount+1;
quiver(x,y,u.',v.',1)
pause(0.01)
 end


xlabel('x')
ylabel('y')
colorbar
title(['stepcount   ' , num2str(stepcount)])
