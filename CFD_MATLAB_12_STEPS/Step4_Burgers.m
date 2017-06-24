clc
clear all
syms x nu t            % Using the symbolic math toolbox
phi = exp(-(x-4*t)^2/(4*nu*(t+1))) + exp(-(x-4*t-2*pi)^2/(4*nu*(t+1)));
phiprime=diff(phi,x); %differentiate with symbolic math enabled
usym = -2*nu*(phiprime/phi)+4; 
ufunc=matlabFunction(usym); %Convert back to Matlab Function
nx = 101;
nt = 100;    %nt is the number of timesteps we want to calculate
nu = 0.07;
dx = 2*pi/(nx-1);
dt = dx*nu;  %dt is the amount of time each timestep covers (delta t)
x=linspace(0,2*pi,nx);
un=zeros(nx);
t=0;



u=ufunc(nu,t,x);


un = zeros(1,nx); %initialize our placeholder array un, to hold the time-stepped solution

for n=1:1:nt;  %iterate through time
    un = u; %%copy the existing values of u into un
    for i=2:1:nx-1  %%now we'll iterate through the u array
    
     %%%This is the line from Step 1, copied exactly.  Edit it for our new equation.
     %%%then uncomment it and run the cell to evaluate Step 2   
      
          % u(i) = un(i)-c*dt/dx*(un(i)-un(i-1)) 
u(i)=un(i)- un(i) * dt/dx *(un(i)-un(i-1))+nu* dt/dx^2 * (un(i+1)-2*un(i) +un(i-1));

    end
u(1) =un(1)- un(1) * dt/dx *(un(1)-un(nx-1))+nu* dt/dx^2 * (un(2)-2*un(1) +un(nx-1)) %Anfangswert 
u(nx)=u(1) %Endwert
    plot(x,u) %%Plot the results
    pause(0.001)
    hold on
    end
    u_analytical = ufunc(nu,n*dt,x);

    plot(x,u_analytical)
    pause(0.001)