% 1D advection test program
% Assumption: Nodal polynomials with node points at interval end points
% Constant transport speed a

clear all

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultTextFontSize',14)

n = 20;         % number of elements
Tf = 5;         % final time
periodic = 0;   % switch between Dirichlet conditions (0) and periodic (1)
k = 1;          % polynomial degree
a = +1;         % advection speed
Cr = 0.5/k^2;   % Courant number -> sets time step size in dt = Cr * h / a
alpha = 0.0;    % flux type, 0 = upwind, 1 = central
left = 0;       % left end of the domain
right = 1;      % right end of the domain
time_integrator = 'rk4'; % rk4, bw_euler, trapezoidal, fw_euler
plot_accurate = 1; % plot only on nodes (0) or with more resolution (1)

% analytical solution
analytical = @(x,t)sin(4*pi*(x-a*t));
% analytical = @(x,t)exp(sin(4*pi*(x-a*t)));
%analytical = @(x,t)(abs(2*mod(x-a*t,1)-1));
%analytical = @(x,t)(mod(x-a*t,1)>0.5);

% set quadrature formula and quadrature nodes for integration
[pg,wg] = get_gauss_quadrature(k+1); 

% set the node points for the Lagrange polynomials, 
xunit = get_gauss_lobatto_quadrature(k+1);

% start up the simulation
% create mesh y and compute mesh size y
kp1 = k+1;
y = zeros(2*n,1);
h = zeros(n,1);
for e=1:n
    y(2*e-1) = left + (right-left)*(e-1)/n;
    y(2*e) = left + (right-left)*e/n;
    h(e) = y(2*e)-y(2*e-1);
end

% create the grid of all node points for interpolation
for e=1:n
    x((kp1*e-k):kp1*e) = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*xunit);
end

% compute time step from Cr number, adjust to hit the desired final time
% exactly
dt = Cr * min(h) / abs(a);
NT = round(Tf/dt);
dt = Tf/NT;

disp(['Number of elements: ' num2str(n) ', minimum mesh size: ' ...
    num2str(min(h)) ', time step size: ' num2str(dt) ])

% flux matrix on a single interface
Fe = 0.5*a*[1 1; -1 -1] + 0.5*abs(a)*(1-alpha)*[1 -1; -1 1];

% mass and advection matrices
[values,derivatives] = evaluate_lagrange_basis(xunit, pg);
Me = values * diag(wg) * values';
Se = values * a * diag(wg) * derivatives';

M = sparse(kp1*n,kp1*n);
Minv = sparse(kp1*n,kp1*n);
S = sparse(kp1*n,kp1*n);
F = sparse(kp1*n,kp1*n);
F(1,1) = Fe(2,2);
F(kp1*n,kp1*n) = Fe(1,1);
for e=1:n
    M((kp1*e-k):kp1*e,(kp1*e-k):kp1*e) = 0.5*h(e)*Me;
    Minv((kp1*e-k):kp1*e,(kp1*e-k):kp1*e) = inv(0.5*h(e)*Me);
    S((kp1*e-k):kp1*e,(kp1*e-k):kp1*e) = Se;
    if e<n
        F(kp1*e:(kp1*e+1),kp1*e:(kp1*e+1)) = Fe;
    end
end
Fbound = zeros(kp1*n,2);
if periodic == 0
    Fbound(1,1) = Fe(2,1);
    Fbound(kp1*n,2) = Fe(1,2);
else
    F(1,kp1*n) = Fe(2,1);
    F(kp1*n,1) = Fe(1,2);
end

% set initial condition
u=zeros(kp1*n,NT+1);
u(:,1) = analytical(x,0);

% run the time loop
for m=1:NT
    
    if strcmp(time_integrator, 'rk4')
        ubound = [analytical(0,(m-1)*dt); analytical(1,(m-1)*dt)];
        k1 = Minv*((S'-F)*u(:,m) - Fbound*ubound);
        ubound = [analytical(0,(m-0.5)*dt); analytical(1,(m-0.5)*dt)];
        k2 = Minv*((S'-F)*(u(:,m)+0.5*dt*k1) - Fbound*ubound);
        k3 = Minv*((S'-F)*(u(:,m)+0.5*dt*k2) - Fbound*ubound);
        ubound = [analytical(0,m*dt); analytical(1,m*dt)];
        k4 = Minv*((S'-F)*(u(:,m)+dt*k3) - Fbound*ubound);
        u(:,m+1) = u(:,m) + dt/6*(k1+2*k2+2*k3+k4);
        
    elseif strcmp(time_integrator, 'fw_euler')
        ubound = [analytical(0,m*dt); analytical(1,m*dt)];
        u(:,m+1) = u(:,m) + dt*Minv*((S'-F)*u(:,m) - Fbound*ubound);
        
    elseif strcmp(time_integrator, 'bw_euler')
        ubound = [analytical(0,m*dt); analytical(1,m*dt)];
        u(:,m+1) = (M/dt - S' + F)\(M/dt * u(:,m) - Fbound*ubound);
        
    elseif strcmp(time_integrator, 'trapezoidal')
        ubound = [analytical(0,m*dt); analytical(1,m*dt)];
        uboundn = [analytical(0,(m-1)*dt); analytical(1,(m-1)*dt)];
        u(:,m+1) = (M/dt - 0.5 * (S'-F))\((M/dt+0.5*(S'-F))*u(:,m) - ...
            0.5*Fbound*ubound - 0.5*Fbound*uboundn);
        
    else
        disp(['unknown integrator ' time_integrator])
        runtime_error;
    end
end

% plot the numerical solution (red), the analytical solution (blue), and
% the initial condition
figure(1)
if plot_accurate == 1
    xx_unit = -1:0.05:1;
    xx=zeros(n,length(xx_unit));
    uu_0 = zeros(size(xx));
    uu = zeros(size(xx));
    val = evaluate_lagrange_basis(xunit, xx_unit);
    for e=1:n
        xx(e,:) = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*xx_unit);
        uu_0(e,:) = val' * u(e*kp1-k:e*kp1,1);
        uu(e,:) = val' * u(e*kp1-k:e*kp1,end);
    end
    u_anal = analytical(xx, Tf);
    % plot dummy data to create correct legend despite multiple columns
    %plot(xx(1,:),uu_0(1,:),'k:',xx(1,:),uu(1,:),'r-',xx(1,:),u_anal(1,:),'b');
    plot(xx(1,:),uu(1,:),'r-',xx(1,:),u_anal(1,:),'b:');
    % plot actual content
    hold on
    %plot(xx',uu_0','k:');
    plot(xx',uu','r-');
    plot(xx',u_anal','b:');
    hold off
else
    plot(x,u(:,1),'k:',x,u(:,end),'r-',x,analytical(x,Tf),'b')
end
xlabel('x')
ylabel('u_h(x)')
title(['degree=' num2str(k) ', n=' num2str(n) ' elements, dt = ' num2str(dt)])
%legend('u_h(x,0)',['u_h(x,' num2str(Tf) ')'],['u(x,' num2str(Tf) ')'])
legend(['u_h(x,' num2str(Tf) ')'],['u(x,' num2str(Tf) ')'])
hold off

% compute error
l2error = 0;
linfty_error = 0;
[pg_err,wg_err] = get_gauss_quadrature(k+3);
values_err = evaluate_lagrange_basis(xunit, pg_err);
for e=1:n
    sol_num = values_err' * u((e-1)*kp1+1:e*kp1,end);
    x_err = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*pg_err);
    sol_exact = analytical(x_err, Tf);
    l2error = l2error + h(e)/2 * wg_err' * (sol_num-sol_exact).^2;
    linfty_error = max([linfty_error; abs(sol_num-sol_exact)]);
end
l2error = sqrt(l2error);
    
disp(['Error in maximum norm ' num2str(linfty_error) ' in L2 norm ' num2str(l2error)])

% plot the solution as it evolves in time
% for n=1:10:N
%     plot(x,v(:,1),'b',x,v(:,n+1),'r-')
%     xlabel('x')
%     ylabel('u_h(x)')
%     legend('u_h(x,0)',['u_h(x,' num2str((n+1)*dt) ')'],'Location','NorthEast')
%     axis([0 1 -1.2 1.2])
%     pause(0.05)
% end

return;

% numerically reproduce the dispersion properties (only equal to the real
% dispersion/dissipation curves for alpha=0)

figure(2)
plot(eig(full(Minv*(S'-F))),'.')

[v,d] = eig(full(Minv*(S'-F)/a));
d = diag(d);
% only take eigenvalues with positive imaginary part
dpos = d(imag(d)>=0);
% sort eigenvalues by their real value
[~,ind] = sort(-(dpos));
figure(3)
lx = 0:1/(length(dpos)-1):1;
% plot dissipation
plot(lx,0*lx,'--',lx,real(dpos(ind))/((k+1)*n),'-')
figure(4)
% plot dispersion
plot(lx,pi*lx,'--',lx,imag(dpos(ind))/((k+1)*n),'-')
