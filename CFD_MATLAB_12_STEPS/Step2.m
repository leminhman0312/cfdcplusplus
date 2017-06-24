nx = 41;
dx = 2./(nx-1);
nt = 40;    %nt is the number of timesteps we want to calculate
dt = .0025;  %dt is the amount of time each timestep covers (delta t)
c=1



u = ones(1,nx)      %as before, we initialize u with every value equal to 1.
u(1,0.5/dx : 1/dx+1)=2  %then set u = 2 between 0.5 and 1 as per our I.C.s

un = ones(1,nx) %initialize our placeholder array un, to hold the time-stepped solutio

for n=1:1:nt;  %iterate through time
    un = u %%copy the existing values of u into un
    for i=2:1:nx  %%now we'll iterate through the u array
    
     %%%This is the line from Step 1, copied exactly.  Edit it for our new equation.
     %%%then uncomment it and run the cell to evaluate Step 2   
      
          % u(i) = un(i)-c*dt/dx*(un(i)-un(i-1)) 
u(i)=un(i)-u(i)*dt/dx*(un(i)-un(i-1))
    
    end
    plot(linspace(0,2,nx),u) %%Plot the results
    pause(0.01)
    end
