function advection_solver_integrate
% 1D advection solver based on pointwise evaluation of fluxes and using
% integrals for the element advection term
% Assumption: Nodal polynomials with node points at interval end points

n = 20;         % number of elements
Tf = 2;         % final time
periodic = 1;   % switch between Dirichlet conditions (0) and periodic (1)
k = 4;          % polynomial degree
a = +1;         % advection speed
Cr = 0.4/k^2;   % Courant number -> sets time step size in dt = Cr * h / a
alpha = 0.0;    % flux type, 0 = upwind, 1 = central
nc = k+1;       % number of quadrature points
left = 0;       % left end of the domain
right = 1;      % right end of the domain
plot_accurate = 1; % plot only on nodes (0) or with more resolution (1)

% analytical solution
%analytical = @(x,t)sin(4*pi*(x-a*t));
analytical = @(x,t)exp(sin(4*pi*(x-a*t)));
%analytical = @(x,t)(abs(2*mod(x-a*t,1)-1));
%analytical = @(x,t)(mod(x-a*t,1)>0.5);

% set quadrature formula and quadrature nodes for integration
[pg,wg] = get_gauss_quadrature(nc); 

% set the node points for the Lagrange polynomials
xunit = get_gauss_lobatto_quadrature(k+1);

% start up the simulation

% create mesh y and compute mesh size h
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

% evaluate reference cell polynomials and mass matrix
[values,derivatives] = evaluate_lagrange_basis(xunit, pg);
Me = values * diag(wg) * values';
Minv = sparse(kp1*n,kp1*n);
for e=1:n
    Minv((kp1*e-k):kp1*e,(kp1*e-k):kp1*e) = inv(0.5*h(e)*Me);
end

% set initial condition
u=zeros(kp1*n,NT+1);
u(:,1) = analytical(x,0);

% run time loop
for m=1:NT
    k1 = Minv*evaluate_adve_rhs(u(:,m), a, ...
                 [analytical(0,(m-1)*dt); analytical(1,(m-1)*dt)],...
                 values, derivatives, wg, alpha, periodic);
    k2 = Minv*evaluate_adve_rhs(u(:,m)+0.5*dt*k1, a, ...
                 [analytical(0,(m-0.5)*dt); analytical(1,(m-0.5)*dt)],...
                 values, derivatives, wg, alpha, periodic);
    k3 = Minv*evaluate_adve_rhs(u(:,m)+0.5*dt*k2, a, ...
                 [analytical(0,(m-0.5)*dt); analytical(1,(m-0.5)*dt)],...
                 values, derivatives, wg, alpha, periodic);
    k4 = Minv*evaluate_adve_rhs(u(:,m)+dt*k3, a, ...
                 [analytical(0,m*dt); analytical(1,m*dt)],...
                 values, derivatives, wg, alpha, periodic);
    
    u(:,m+1) = u(:,m) + dt/6*(k1+2*k2+2*k3+k4);
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
    plot(xx(1,:),uu_0(1,:),'k:',xx(1,:),uu(1,:),'r-',xx(1,:),u_anal(1,:),'b');
    % plot actual content
    hold on
    plot(xx',uu_0','k:');
    plot(xx',uu','r-');
    plot(xx',u_anal','b-');
    hold off
else
    plot(x,u(:,1),'k:',x,u(:,end),'r-',x,analytical(x,Tf),'b')
end
xlabel('x')
ylabel('u_h(x)')
title(['degree=' num2str(k) ', n=' num2str(n) ' elements, dt = ' num2str(dt)])
legend('u_h(x,0)',['u_h(x,' num2str(Tf) ')'],['u(x,' num2str(Tf) ')'])

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

end


function rhs = evaluate_adve_rhs(u, a, bc, values, derivatives, weights, alpha, periodic)

kp1 = size(values, 1); % degree + 1
n = length(u)/kp1;
rhs = zeros(size(u));

for e=1:n
    ue = u((e-1)*kp1+1:e*kp1);
    % interpolate u to quadrature points
    u_quad = (values')*ue;
    % compute operator at quadrature points and multiply by gradient of
    % test function
    rhs((e-1)*kp1+1:e*kp1) = derivatives * (a * weights .* u_quad);
    
    % compute advective numerical flux on the left
    uminus = ue(1);
    if (e==1)
        if (periodic)
            % periodic bc
            uplus = u(end);
        else
            % Dirichlet condition, implemented via mirror principle
            % uplus = 2*bc(1)-uminus;
            % Dirichlet condition, implemented in classical way
            uplus = bc(1);
        end
    else
        uplus = u((e-1)*kp1);
    end
    numflux = (uminus*a+uplus*a)/2 + (1-alpha)/2 * abs(a)*(uplus-uminus);
    rhs((e-1)*kp1+1) = rhs((e-1)*kp1+1) + numflux;
    
    % compute advective numerical flux on the right
    uminus = ue(kp1);
    if (e==n)
        if (periodic)
            % periodic bc
            uplus = u(1);
        else
            % Dirichlet condition, implemented via mirror principle
            % uplus = 2*bc(2)-uminus;
            % Dirichlet condition, implemented in classical way
            uplus = bc(2);
        end
    else
        uplus = u(e*kp1+1);
    end
    numflux = (uminus*a+uplus*a)/2 + (1-alpha)/2 * abs(a)*(uminus-uplus);
    rhs(e*kp1) = rhs(e*kp1) - numflux;    
end

end
