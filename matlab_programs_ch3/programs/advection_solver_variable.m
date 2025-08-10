function advection_solver_variable
% 1D advection solver based on pointwise evaluation of fluxes and using
% integrals for the element advection term, using a variable transport
% speed and periodic boundary conditions
% Assumption: Nodal polynomials with node points at interval end points

n = 10;         % number of elements
Tf = 2;         % final time
k = 5;          % polynomial degree
periodic = 1;   % switch between Dirichlet conditions (0) and periodic (1)
Cr = 0.4/k^2;   % Courant number -> sets time step size in dt = Cr * h / a
alpha = 0.0;    % flux type, 0 = upwind, 1 = central
nc = k+1;       % number of quadrature points
left = -1;      % left end of the domain
right = 1;      % right end of the domain
plot_accurate = 1; % plot only on nodes (0) or with more resolution (1)

% variable advection speed
a = @(x,t)1;%(1-x.^2).^5+1;

% solution only valid for a=1
%analytical = @(x,t)sin(4*pi*(x-t));
analytical = @(x,t)exp(sin(10*pi*(x-t)));
%analytical = @(x,t)(abs(2*mod(x-t,1)-1));
%analytical = @(x,t)(mod(x-t,1)>0.5);

% set quadrature formula and quadrature nodes for integration
[pg,wg] = get_gauss_quadrature(nc);

% set the node points for the Lagrange polynomials
xunit = get_gauss_lobatto_quadrature(k+1);

% start up the simulation

% create mesh y and compute mesh size h
kp1 = k+1;
y = zeros(n+1,1);
h = zeros(n,1);
y(1) = left;
for e=1:n
    y(e+1) = left + (right-left)*e/n;
    h(e) = y(e+1)-y(e);
end

% create the grid of all node points for interpolation
for e=1:n
    x((kp1*e-k):kp1*e) = y(e)+(y(e+1)-y(e))*(0.5+0.5*xunit);
end

% compute time step from Cr number, adjust to hit the desired final time
% exactly
dt = Cr * min(h) / max(abs(a(x,0)));
NT = round(Tf/dt);
dt = Tf/NT;

disp(['Number of elements: ' num2str(n) ', minimum mesh size: ' ...
    num2str(min(h)) ', time step size: ' num2str(dt) ])

% evaluate reference cell polynomials and mass matrix
[values,derivatives] = evaluate_lagrange_basis(xunit, pg);
Me = values * diag(wg) * values';
Minv = sparse(kp1*n,kp1*n);

tic

for e=1:n
    Minv((kp1*e-k):kp1*e,(kp1*e-k):kp1*e) = inv(0.5*h(e)*Me);
end

% set initial condition
u=zeros(kp1*n,NT+1);
u(:,1) = analytical(x,0);

% run time loop
for m=1:NT
    k1 = Minv*evaluate_adve_rhs(u(:,m), a, y, (m-1)*dt, ...
                 analytical, values, derivatives, pg, wg, alpha, periodic);
    k2 = Minv*evaluate_adve_rhs(u(:,m)+0.5*dt*k1, a, y, (m-0.5)*dt, ...
                 analytical, values, derivatives, pg, wg, alpha, periodic);
    k3 = Minv*evaluate_adve_rhs(u(:,m)+0.5*dt*k2, a, y, (m-0.5)*dt, ...
                 analytical, values, derivatives, pg, wg, alpha, periodic);
    k4 = Minv*evaluate_adve_rhs(u(:,m)+dt*k3, a, y, m*dt, ...
                 analytical, values, derivatives, pg, wg, alpha, periodic);
    
    u(:,m+1) = u(:,m) + dt/6*(k1+2*k2+2*k3+k4);
end

toc

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
        xx(e,:) = y(e)+(y(e+1)-y(e))*(0.5+0.5*xx_unit);
        uu_0(e,:) = val' * u(e*kp1-k:e*kp1,1);
        uu(e,:) = val' * u(e*kp1-k:e*kp1,end);
    end
    % plot dummy data to create correct legend despite multiple columns
    plot(xx(1,:),uu_0(1,:),'k:',xx(1,:),uu(1,:),'r-');
    % plot actual content
    hold on
    plot(xx',uu_0','k:');
    plot(xx',uu','r-');
    hold off
else
    plot(x,u(:,1),'k:',x,u(:,end),'r-')
end
xlabel('x')
ylabel('u_h(x)')
title(['degree=' num2str(k) ', n=' num2str(n) ' elements, dt = ' num2str(dt)])
legend('u_h(x,0)',['u_h(x,' num2str(Tf) ')'])

l2error = 0;
linfty_error = 0;
[pg_err,wg_err] = get_gauss_quadrature(k+3);
values_err = evaluate_lagrange_basis(xunit, pg_err);
mass_init = 0;
mass_end = 0;
energy_init = 0;
energy_end = 0;
for e=1:n
    sol_num = values_err' * u((e-1)*kp1+1:e*kp1,end);
    
    mass_init = mass_init + h(2)/2 * wg_err' * ...
        (values_err' * u((e-1)*kp1+1:e*kp1,1));
    mass_end = mass_end + h(e)/2 * wg_err' * sol_num;
    
    energy_init = energy_init + h(2)/2 * wg_err' * ...
        (values_err' * u((e-1)*kp1+1:e*kp1,1)).^2;
    energy_end = energy_end + h(e)/2 * wg_err' * sol_num.^2;
    
    x_err = y(e)+(y(e+1)-y(e))*(0.5+0.5*pg_err);
    sol_exact = analytical(x_err, Tf);
    l2error = l2error + h(e)/2 * wg_err' * (sol_num-sol_exact).^2;
    linfty_error = max([linfty_error; abs(sol_num-sol_exact)]);
end
l2error = sqrt(l2error);

disp(['Distance in maximum norm ' num2str(linfty_error) ' in L2 norm ' num2str(l2error)])
disp(['Mass t=0: ' num2str(mass_init) ' mass at t=' num2str(Tf) ': '...
    num2str(mass_end) ' difference ' num2str(mass_init-mass_end)])
disp(['Energy t=0: ' num2str(energy_init) ' energy at t=' num2str(Tf) ': '...
    num2str(energy_end) ' difference ' num2str(energy_init-energy_end)])

end


function rhs = evaluate_adve_rhs(u, a, y, t, analytical, values, ...
    derivatives, pg, weights, alpha, periodic)

kp1 = size(values, 1); % degree + 1
n = length(u)/kp1;
rhs = zeros(size(u));
normal = [-1 1];
face_indices = [1 kp1];
neighbor_index = [-1 1];
neighbor_face_indices = [kp1 1];

for e=1:n
    ue = u((e-1)*kp1+1:e*kp1);
    xe = y(e) + (y(e+1)-y(e))*(0.5 + 0.5 * pg);
    ae = a(xe,t);
    % interpolate u to quadrature points
    u_quad = (values')*ue;
    % compute operator at quadrature points and multiply by gradient of
    % test function
    myrhs = derivatives * (ae .* weights .* u_quad);
    
    % compute advective numerical flux on the two faces
    for f=1:2
        uminus = ue(face_indices(1,f));
        aminus = a(y(e+f-1)-normal(f)*1e-13,t);
        if (f==1 && e==1) || (f==2 && e==n)
            if periodic
                uplus = u((2-f)*(n-1)*kp1+neighbor_face_indices(f));
            else
                uplus = analytical(y(1+(f-1)*n),t);
            end
            aplus = aminus;
        else
            uplus = u((e-1+neighbor_index(f))*kp1+neighbor_face_indices(f));
            aplus = a(y(e+f-1)+normal(f)*1e-13,t);
        end
        numflux = (uminus*aminus+uplus*aplus)/2*normal(f) + ...
                  (1-alpha)/2 * max(abs(aminus),abs(aplus))*(uminus-uplus);
        myrhs(face_indices(1,f)) = myrhs(face_indices(1,f)) - numflux;
    end
    rhs((e-1)*kp1+1:e*kp1) = myrhs;
end

end
