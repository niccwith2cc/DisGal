% 1D advection test program implementing a 2nd order central finite
% difference method for constant transport speed a

%clear all

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultTextFontSize',14)

n = 321;         % number of grid points
Tf = 5;         % final time
periodic = 0;   % switch between Dirichlet conditions (0) and periodic (1)
a = +1;         % advection speed
Cr = 0.4;       % Courant number -> sets time step size in dt = Cr * h / a
left = 0;       % left end of the domain
right = 1;      % right end of the domain
time_integrator = 'rk4'; % rk4, bw_euler, trapezoidal, fw_euler

% analytical solution
analytical = @(x,t)sin(4*pi*(x-a*t));
%analytical = @(x,t)exp(sin(4*pi*(x-a*t)));

h = (right-left)/(n-1);
% this is the basic central difference formula in the interior
D = -a / 2 / h * spdiags([-ones(n,1) zeros(n,1) ones(n,1)], -1:1, n, n);
% one-sided finite difference formula at the boundary & simultaneous
% approximation term for the boundary condition, i.e., a penalty 
% tau * (u_1 - h_0) for tau = 2/h
D(1,2) = -a / h;
D(1,1) = a / h - (a + abs(a)) / h;
D(n,n-1) = a / h;
D(n,n) = -a / h + (a - abs(a)) / h;

x = left:h:right;

Fbound = zeros(n,2);
if periodic == 0
    Fbound(1,1) = -(a + abs(a)) / h;
    Fbound(n,2) =  (a - abs(a)) / h;
else
    D(1,n) =  (a + abs(a)) / h;
    D(n,1) = -(a - abs(a)) / h;
end

% compute time step from Cr number, adjust to hit the desired final time
% exactly
dt = Cr * min(h) / abs(a);
NT = round(Tf/dt);
dt = Tf/NT;

% set initial condition
u=zeros(n,NT+1);
u(:,1) = analytical(x,0);

% run the time loop
for m=1:NT
    
    if strcmp(time_integrator, 'rk4')
        ubound = [analytical(0,(m-1)*dt); analytical(1,(m-1)*dt)];
        k1 = D*u(:,m) - Fbound*ubound;
        ubound = [analytical(0,(m-0.5)*dt); analytical(1,(m-0.5)*dt)];
        k2 = D*(u(:,m)+0.5*dt*k1) - Fbound*ubound;
        k3 = D*(u(:,m)+0.5*dt*k2) - Fbound*ubound;
        ubound = [analytical(0,m*dt); analytical(1,m*dt)];
        k4 = D*(u(:,m)+dt*k3) - Fbound*ubound;
        u(:,m+1) = u(:,m) + dt/6*(k1+2*k2+2*k3+k4);
        
    elseif strcmp(time_integrator, 'fw_euler')
        ubound = [analytical(0,m*dt); analytical(1,m*dt)];
        u(:,m+1) = u(:,m) + dt*(D*u(:,m) - Fbound*ubound);
        
    elseif strcmp(time_integrator, 'bw_euler')
        ubound = [analytical(0,m*dt); analytical(1,m*dt)];
        u(:,m+1) = (eye(n)/dt - D)\(1/dt * u(:,m) - Fbound*ubound);
        
    elseif strcmp(time_integrator, 'trapezoidal')
        ubound = [analytical(0,m*dt); analytical(1,m*dt)];
        uboundn = [analytical(0,(m-1)*dt); analytical(1,(m-1)*dt)];
        u(:,m+1) = (eye(n)/dt - 0.5 * D)\((eye(n)/dt+0.5*D)*u(:,m) - ...
            0.5*Fbound*ubound - 0.5*Fbound*uboundn);
        
    else
        disp(['unknown integrator ' time_integrator])
        runtime_error;
    end
end

% plot the initial condition and the final solution
figure(1)
plot(x,u(:,1),'b',x,u(:,end),'r-')
xlabel('x')
ylabel('u_h(x)')
title(['n=' num2str(n) ', dt = ' num2str(dt)])
legend('u_h(x,0)',['u_h(x,' num2str(Tf) ')'],'Location','NorthEast')

error = max(abs(u(:,end) - analytical(x',Tf)));
disp(['Error in maximum norm: ' num2str(error)]);