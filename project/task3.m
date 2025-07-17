% function acoustic_solver_integrate
% 1D acoustic solver based on pointwise evaluation of fluxes and using
% integrals for the element advection term
% Assumption: Nodal polynomials with node points at interval end points
close all
clear

n = 20;         % number of elements
Tf = 3e-3;      % final time
bc = 0;         % select boundary condition: Dirichlet (0), Neumann(1), Absorbing(2)
k = 2;          % polynomial degree
c = 340;        % advection speed
rho = 1.2;      % density
Cr = 0.4/k^2;   % Courant number -> sets time step size in dt = Cr * h / a
alpha = 0.0;    % flux type, 0 = upwind, 1 = central
nc = k+1;       % number of quadrature points
left = 0;       % left end of the domain
right = 1;      % right end of the domain
plot_accurate = 1; % plot only on nodes (0) or with more resolution (1)
flux_type = 1;  % choose between Lax-Friedrich (0) or HDG (1)

% analytical solution
analytical_v = @(x,t)zeros(size(x));
analytical_p = @(x,t)exp(-((x - 0.5).^2) / 0.02^2);

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
x = zeros(1, kp1*n);
for e=1:n
    x((e-1)*kp1+1:e*kp1) = y(2*e-1) + h(e)/2*(xunit + 1);
end

% compute time step from Cr number, adjust to hit the desired final time
% exactly
dt = Cr * min(h) / abs(c);
NT = round(Tf/dt);
dt = Tf/NT;

disp(['Number of elements: ' num2str(n) ', polynomial degree: ' num2str(k) ', minimum mesh size: ' ...
    num2str(min(h)) ', time step size: ' num2str(dt) ])

switch flux_type
    case 0
        flux_name = 'Lax';
    case 1
        flux_name = 'HDG';
end
switch bc
    case 0
        bc_name = 'Dirichlet';
    case 1
        bc_name = 'Neumann';
    case 2
        bc_name = 'Absorbing';
end
disp(['Flux used: ', flux_name, ' with ', bc_name, ' BC']);

tic;
% evaluate reference cell polynomials and mass matrix
[values,derivatives] = evaluate_lagrange_basis(xunit, pg);
Me = values * diag(wg) * values';
Minv = sparse(2*kp1*n, 2*kp1*n);
for e=1:n
    Mloc = inv(0.5 * h(e) * Me);
    Minv( (kp1*e-k):kp1*e, (kp1*e-k):kp1*e ) = Mloc; % for v
    Minv( (e+n)*kp1-k:(e+n)*kp1, (e+n)*kp1-k:(e+n)*kp1 ) = Mloc; % for p
end

% set initial conditions
w = zeros(2*kp1*n, NT+1);
w0 = [analytical_v(x, 0), analytical_p(x, 0)];
w(:, 1) = w0;

% run time loop
for m=1:NT
    k1 = Minv * evaluate_acoustic_rhs(flux_type, w(:, m), c, rho, [analytical_v(0, (m-1)*dt); analytical_v(1, (m-1)*dt); analytical_p(0, (m-1)*dt); analytical_p(1, (m-1)*dt)], values, derivatives, wg, bc);

    k2 = Minv * evaluate_acoustic_rhs(flux_type, w(:, m) + 0.5*dt*k1, c, rho, [analytical_v(0, (m-0.5)*dt); analytical_v(1, (m-0.5)*dt); analytical_p(0, (m-0.5)*dt); analytical_p(1, (m-0.5)*dt)], values, derivatives, wg, bc);
                 
    k3 = Minv * evaluate_acoustic_rhs(flux_type, w(:, m) + 0.5*dt*k2, c, rho, [analytical_v(0, (m-0.5)*dt) ; analytical_v(1, (m-0.5)*dt); analytical_p(0, (m-0.5)*dt) ; analytical_p(1, (m-0.5)*dt)], values, derivatives, wg, bc);

    k4 = Minv * evaluate_acoustic_rhs(flux_type, w(:, m) + dt*k3, c, rho, [analytical_v(0, m*dt); analytical_v(1, m*dt); analytical_p(0, m*dt); analytical_p(1, m*dt)], values, derivatives, wg, bc);
    
    w(:, m+1) = w(:, m) + dt/6*(k1+2*k2+2*k3+k4);
end

% calculate L2 and Linf for velocity
l2error_v = 0;
linfty_error_v = 0;
[pg_err,wg_err] = get_gauss_quadrature(k+3);
values_err = evaluate_lagrange_basis(xunit, pg_err);
for e=1:n
    sol_num = values_err' * w((e-1)*kp1+1:(e)*kp1, end);
    x_err = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*pg_err);
    sol_exact = analytical_v(x_err, Tf);
    l2error_v = l2error_v + h(e)/2 * wg_err' * (sol_num-sol_exact).^2;
    linfty_error_v = max([linfty_error_v; abs(sol_num-sol_exact)]);
end
l2error_v = sqrt(l2error_v);

disp(['Velocity Error in maximum norm ' num2str(linfty_error_v) ' in L2 norm ' num2str(l2error_v)])

% calculate L2 and Linf for pressure
l2error_p = 0;
linfty_error_p = 0;
[pg_err,wg_err] = get_gauss_quadrature(k+3);
values_err = evaluate_lagrange_basis(xunit, pg_err);
for e=1:n
    sol_num = values_err' * w((e+n-1)*kp1+1:(e+n)*kp1, end);
    x_err = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*pg_err);
    sol_exact = analytical_p(x_err, Tf);
    l2error_p = l2error_p + h(e)/2 * wg_err' * (sol_num-sol_exact).^2;
    linfty_error_p = max([linfty_error_p; abs(sol_num-sol_exact)]);
end
l2error_p = sqrt(l2error_p);

disp(['Pressure Error in maximum norm ' num2str(linfty_error_p) ' in L2 norm ' num2str(l2error_p)])

toc;

% plotting solution from 0 to tf every 6e-4
dt_plot = 6e-4;
NT_plot = round(Tf/dt_plot);
t_plot = 0:dt_plot:Tf;
index_t = round((t_plot/Tf)*NT) + 1;

figure;
xx_unit = -1:0.05:1;
xx = zeros(n, length(xx_unit));
vvv = zeros(size(xx, 1), size(xx, 2), NT_plot+1); % creates a 3D matrix of size(xx), size(xx), NT_plot+1
val = evaluate_lagrange_basis(xunit, xx_unit);
for e=1:n
    xx(e, :) = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*xx_unit);
    vvv(e, :, 1) = val' * w((e)*kp1-k: (e)*kp1, index_t(1));
    vvv(e, :, 2) = val' * w((e)*kp1-k: (e)*kp1, index_t(2));
    vvv(e, :, 3) = val' * w((e)*kp1-k: (e)*kp1, index_t(3));
    vvv(e, :, 4) = val' * w((e)*kp1-k: (e)*kp1, index_t(4));
    vvv(e, :, 5) = val' * w((e)*kp1-k: (e)*kp1, index_t(5));
    vvv(e, :, 6) = val' * w((e)*kp1-k: (e)*kp1, index_t(6));
end

v_anal = zeros(size(xx, 1), size(xx, 2), NT_plot+1);
for i = 1:(NT_plot+1)
    v_anal(:,:,i) = analytical_v(xx, t_plot(i));
end

for i = 1:6
    subplot(2, 3, i);
    plot(xx(1,:), vvv(1, :, i),'r-', xx(1,:), v_anal(1, :, i),'b');
    title(['Subplot ', num2str(i)]);
    hold on
    plot(xx',vvv(:,:,i)','r-', 'LineWidth', 1.3);
    plot(xx',v_anal(:,:,i)','b-');
    hold off
    xlabel('x')
    ylabel('v_h(x)')
    title(['t = ' num2str(t_plot(i))])
    legend(['v_h(x, ' num2str(t_plot(i)) ')'], 'v(x, 0)')
end

sgtitle(['Velocity evolution with ' flux_name ', BC = ' bc_name ', degree = ' num2str(k) ', elements = ' num2str(n)]);

figure;
xx_unit = -1:0.05:1;
xx = zeros(n, length(xx_unit));
ppp = zeros(size(xx, 1), size(xx, 2), NT_plot+1); % creates a 3D matrix of size(xx), size(xx), NT_plot+1
val = evaluate_lagrange_basis(xunit, xx_unit);
for e=1:n
    xx(e, :) = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*xx_unit);
    ppp(e, :, 1) = val' * w((e+n)*kp1-k: (e+n)*kp1, index_t(1));
    ppp(e, :, 2) = val' * w((e+n)*kp1-k: (e+n)*kp1, index_t(2));
    ppp(e, :, 3) = val' * w((e+n)*kp1-k: (e+n)*kp1, index_t(3));
    ppp(e, :, 4) = val' * w((e+n)*kp1-k: (e+n)*kp1, index_t(4));
    ppp(e, :, 5) = val' * w((e+n)*kp1-k: (e+n)*kp1, index_t(5));
    ppp(e, :, 6) = val' * w((e+n)*kp1-k: (e+n)*kp1, index_t(6));
end

p_anal = zeros(size(xx, 1), size(xx, 2), NT_plot+1);
for i = 1:(NT_plot+1)
    p_anal(:,:,i) = analytical_p(xx, t_plot(i));
end

for i = 1:6
    subplot(2, 3, i);
    plot(xx(1,:), ppp(1, :, i),'r-', xx(1,:), p_anal(1, :, i),'b');
    title(['Subplot ', num2str(i)]);
    hold on
    plot(xx',ppp(:,:,i)','r-', 'LineWidth', 1.3);
    plot(xx',p_anal(:,:,i)','b-');
    hold off
    xlabel('x')
    ylabel('p_h(x)')
    title(['t = ' num2str(t_plot(i))])
    legend(['p_h(x, ' num2str(t_plot(i)) ')'], 'p(x, 0)')
end

sgtitle(['Pressure evolution with ' flux_name ', BC = ' bc_name ', degree = ' num2str(k) ', elements = ' num2str(n)]);

function rhs = evaluate_acoustic_rhs(flux_type, w, c, rho, bv, values, derivatives, weights, bc)
    switch flux_type
        case 0 % local Lax-Friedrich flux
            rhs = lax_flux_rhs(w, c, rho, bv, values, derivatives, weights, bc);
        
        case 1 % HDG flux
            rhs = hdg_flux_rhs(w, c, rho, bv, values, derivatives, weights, bc);
    end
end

function rhs = lax_flux_rhs(w, c, rho, bv, values, derivatives, weights, bc)
% ------------------------------------------------------------------------------------------- %
% function to calculate rhs using local Lax-Friedrich flux
% VARIABLES:
%   IN:
% w           : system of equation of w = [v, p] ---- size of 2*n*kp1
% c           : speed of sound
% rho         : density of fluid
% bv          : boundary values with respect to the boundary conditions in the form of [vL; vR; pL; pR]
% values      : lagrange polynomial values
% derivatives : lagrange polynomial derivatives
% weights     : weights of the gauss quadrature
% bc          : Boundary condition
% ------------------------------------------------------------------------------------------- %
%   OUT:
% rhs: vector of right hand side to be multiplied by Minv
% ------------------------------------------------------------------------------------------- %

kp1 = size(values, 1); % degree + 1
n = length(w)/(2*kp1);
rhs = zeros(size(w));

for e=1:n
    ve = w((e-1)*kp1+1: e*kp1); % 1 to kp1*n
    pe = w((e+n-1)*kp1+1: (e+n)*kp1); % from kp1*n+1 to 2*kp1*n

    % interpolate v to quadrature points
    v_quad = values' * ve;
    p_quad = values' * pe;
    flux_v = weights .* p_quad;
    flux_p = weights .* v_quad;
    
    % compute operator at quadrature points and multiply by gradient of
    % test function
    rhs((e-1)*kp1+1: e*kp1) = derivatives * flux_v / rho;
    rhs((e+n-1)*kp1+1: (e+n)*kp1) = derivatives * flux_p * rho * c^2;
    
    % compute acoustic numerical flux on the left
    vminus = ve(1);
    pminus = pe(1);
    if (e==1)
        switch bc
            case 0 % Dirichlet condition
                vplus = vminus;
                pplus = 2*bv(3)-pminus;
            case 1 % Neumann
                vplus = -bv(1)/rho;
                pplus = pminus;
            case 2 % Absorbing
                vplus = -pminus/(rho*c);
                pplus = pminus;
        end
    else
        vplus = w((e-1)*kp1); % e > 1: w(kp1:(e-1)*kp1)
        pplus = w((e+n-1)*kp1); 
    end
    numflux_v = 1/(2*rho) * (pminus + pplus) + 1/2 * abs(c) * (vplus - vminus);
    numflux_p = 1/2*c^2*rho * (vminus + vplus) + 1/2 * abs(c) * (pplus - pminus);
    rhs((e-1)*kp1+1) = rhs((e-1)*kp1+1) + numflux_v;
    rhs((e+n-1)*kp1+1) = rhs((e+n-1)*kp1+1) + numflux_p;
    
    % compute acoustic numerical flux on the right
    vminus = ve(kp1);
    pminus = pe(kp1);
    if (e==n)
        switch bc
            case 0 % Dirichlet condition
                vplus = vminus;
                pplus = 2*bv(4)-pminus;
            case 1 % Neumann
                vplus = bv(2)/rho;
                pplus = pminus;
            case 2 % Absorbing
                vplus = pminus/(rho*c);
                pplus = pminus;
        end
    else
        vplus = w(e*kp1+1);
        pplus = w((e+n)*kp1+1);
    end
    numflux_v  = 1/(2*rho) * (pminus + pplus) + 1/2 * abs(c)*(vminus - vplus);
    numflux_p = 1/2*c^2*rho * (vminus + vplus) + 1/2 * abs(c)*(pminus - pplus);
    rhs(e*kp1) = rhs(e*kp1) - numflux_v;
    rhs((e+n)*kp1) = rhs((e+n)*kp1) - numflux_p;
end

end

function rhs = hdg_flux_rhs(w, c, rho, bv, values, derivatives, weights, bc)
% ------------------------------------------------------------------------------------------- %
% function to calculate rhs using HDG flux
% VARIABLES:
%   IN:
% w           : system of equation of w = [v, p] ---- size of 2*n*kp1
% c           : speed of sound
% rho         : density of fluid
% bv          : boundary values with respect to the boundary conditions in the form of [vL; vR; pL; pR]
% values      : lagrange polynomial values
% derivatives : lagrange polynomial derivatives
% weights     : weights of the gauss quadrature
% bc          : Boundary condition
% ------------------------------------------------------------------------------------------- %
%   OUT:
% rhs: vector of right hand side to be multiplied by Minv
% ------------------------------------------------------------------------------------------- %

kp1 = size(values, 1); % degree + 1
n = length(w)/(2*kp1);
rhs = zeros(size(w));

% stabilization parameters
tau = 1/c;

for e = 1:n
    ve = w((e-1)*kp1+1:e*kp1);
    pe = w((e+n-1)*kp1+1:(e+n)*kp1);

    % Evaluate at quadrature
    v_quad = values' * ve;
    p_quad = values' * pe;
    flux_v = weights .* p_quad;
    flux_p = weights .* v_quad;

    rhs((e-1)*kp1+1:e*kp1) = derivatives * flux_v/rho;
    rhs((e+n-1)*kp1+1:(e+n)*kp1) = derivatives * flux_p*rho*c^2;

    % == LEFT INTERFACE ==
    vminus = ve(1);
    pminus = pe(1);
    norm = -1;
    if e == 1
        switch bc
            case 0 % Dirichlet
                pplus = bv(3);
                lambda = pplus;
            case 1 % Neumann check n
                lambda = (1/tau) * rho * vminus * norm + pminus;
            case 2 % Absorbing check n
                lambda = (1 /(tau + 1/c))* (rho * vminus * norm + tau * pminus);
        end
    else 
        vplus = w((e-1)*kp1); % e > 1: w(kp1:(e-1)*kp1)
        pplus = w((e+n-1)*kp1);
        % check n
        lambda = (rho / (2*tau)) * (vminus - vplus) * norm + 1/2 * (pminus + pplus);
    end

    numflux_v = lambda/rho;
    % check n
    numflux_p = (vminus + (tau/rho) * norm * (pminus - lambda))*rho*c^2;

    rhs((e-1)*kp1+1) = rhs((e-1)*kp1+1) + numflux_v;
    rhs((e+n-1)*kp1+1) = rhs((e+n-1)*kp1+1) + numflux_p;

    % == RIGHT INTERFACE ==
    vminus = ve(kp1); 
    pminus = pe(kp1);
    norm = +1;
    if e == n
        switch bc
            case 0 % Dirichlet
                pplus = bv(4);
                lambda = pplus;
            case 1 % Neumann check n
                lambda = (1/tau) * rho * vminus * norm + pminus;
            case 2 % Absorbing check n
                lambda = (1 /(tau + 1/c)) * (rho * vminus * norm + tau * pminus);
        end
    else
        vplus = w(e*kp1+1); % e > 1: w(kp1:(e-1)*kp1)
        pplus = w((e+n)*kp1+1); 
        % check n
        lambda = (rho / (2*tau)) * (vminus - vplus) * norm + 1/2 * (pminus + pplus);
    end

    numflux_v = lambda/rho;
    % check n
    numflux_p = (vminus + (tau/rho) * norm * (pminus - lambda))*rho*c^2;
    rhs(e*kp1) = rhs(e*kp1) - numflux_v;
    rhs((e+n)*kp1) = rhs((e+n)*kp1) - numflux_p;
end
end