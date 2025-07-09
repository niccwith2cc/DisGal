% 1D advection test program for continuous finite elements with boundary
% conditions implemented strongly. Only implicit time integration with
% the backward Euler method and the trapezoidal rule are implemented.
% The program assumes a constant transport speed a > 0

%clear all

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultTextFontSize',14)

n = 40;         % number of elements
Tf = 5;         % final time
periodic = 0;   % switch between Dirichlet conditions (0) and periodic (1)
k = 1;          % polynomial degree
a = +1;         % advection speed, only positive supported
Cr = 0.2/k;     % Courant number -> sets time step size in dt = Cr * h / a
left = 0;       % left end of the domain
right = 1;      % right end of the domain
time_integrator = 'trapezoidal'; % bw_euler, trapezoidal

% analytical solution
analytical = @(x,t)sin(4*pi*(x-a*t));
%analytical = @(x,t)exp(sin(4*pi*(x-a*t)));

% set quadrature formula and quadrature nodes for integration
[pg,wg] = get_gauss_quadrature(k+1); 

% set the node points for the Lagrange polynomials, 
xunit = get_gauss_lobatto_quadrature(k+1);

% start up the simulation
% create mesh y and compute mesh size y
y = zeros(2*n,1);
h = zeros(n,1);
for e=1:n
    y(2*e-1) = left + (right-left)*(e-1)/n;
    y(2*e) = left + (right-left)*e/n;
    h(e) = y(2*e)-y(2*e-1);
end

% create the grid of all node points for interpolation
for e=1:n
    x((k*e-k+1):k*e+1) = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*xunit);
end

% compute time step from Cr number, adjust to hit the desired final time
% exactly
dt = Cr * min(h) / abs(a);
NT = round(Tf/dt);
dt = Tf/NT;

disp(['Number of elements: ' num2str(n) ', minimum mesh size: ' ...
    num2str(min(h)) ', time step size: ' num2str(dt) ])

% mass and advection matrices
[values,derivatives] = evaluate_lagrange_basis(xunit, pg);
Me = values * diag(wg) * values';
Se = values * a * diag(wg) * derivatives';

M = sparse(k*n+1,k*n+1);
S = sparse(k*n+1,k*n+1);
for e=1:n
    M((k*e-k+1):k*e+1,(k*e-k+1):k*e+1) = ...
        M((k*e-k+1):k*e+1,(k*e-k+1):k*e+1) + 0.5*h(e)*Me;
    S((k*e-k+1):k*e+1,(k*e-k+1):k*e+1) = ...
        S((k*e-k+1):k*e+1,(k*e-k+1):k*e+1) + Se;
end
if periodic == 0
    % Dirichlet: set first entry in matrix to zero
    M(1,1) = 1;
    M(1,2:k+1) = 0;
    S(1,1) = 1;
    S(1,2:k+1) = 0;
else
    % periodic: eliminate last degree of freedom in favor of first one
    M(1,k*n-k+1:k*n) = M(k*n+1,k*n-k+1:k*n);
    M(k*n-k+1:k*n,1) = M(k*n-k+1:k*n,k*n+1);
    M(1,1) = M(1,1) + M(end,end);
    M = M(1:end-1,1:end-1);
    S(1,k*n-k+1:k*n) = S(k*n+1,k*n-k+1:k*n);
    S(k*n-k+1:k*n,1) = S(k*n-k+1:k*n,k*n+1);
    S(1,1) = S(1,1) + S(end,end);
    S = S(1:end-1,1:end-1);
    x = x(1:end-1);
end

% set initial condition
u=zeros(size(M,1),NT+1);
u(:,1) = analytical(x,0);

% run the time loop
for m=1:NT
    
    if strcmp(time_integrator, 'bw_euler')
        rhs = M/dt * u(:,m);
        if periodic == 0
            rhs(1) = (M(1,1)/dt + S(1,1)) * analytical(0,m*dt);
        end
        u(:,m+1) = (M/dt + S)\rhs;
        
    elseif strcmp(time_integrator, 'trapezoidal')
        rhs = (M/dt - 0.5*S)*u(:,m);
        if periodic == 0
            rhs(1) = (M(1,1)/dt + 0.5*S(1,1)) * analytical(0,m*dt);
        end
        u(:,m+1) = (M/dt + 0.5 * S)\rhs;
        
    else
        disp(['unknown integrator ' time_integrator])
        runtime_error;
    end
end

% plot the initial condition and the final solution, fill up in case of
% periodic boundary conditions
if periodic == 1
    x = [x y(end)];
    u = [u; u(1,:)];
end
figure(1)
plot(x,analytical(x,Tf),'b',x,u(:,end),'r-')
xlabel('x')
ylabel('u_h(x)')
title(['degree=' num2str(k) ', n_{ele}=' num2str(n) ', dt = ' num2str(dt)])
legend('u_h(x,0)',['u_h(x,' num2str(Tf) ')'],'Location','NorthEast')

% compute error
l2error = 0;
linfty_error = 0;
[pg_err,wg_err] = get_gauss_quadrature(k+3);
values_err = evaluate_lagrange_basis(xunit, pg_err);
for e=1:n
    sol_num = values_err' * u((e-1)*k+1:e*k+1,end);
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

return