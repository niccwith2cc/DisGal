close all
clear
Courant = 1:0.01:3.08;
quad_type = [0, 1]; % quadrature type (0) Gauss, (1) Gauss Lobatto

soft_purple = [0.6, 0.4, 0.8];
dusty_pink = [0.9, 0.6, 0.7];
muted_green = [0.4, 0.7, 0.6];
warm_red = [0.8, 0.4, 0.4];
soft_blue = [0.4, 0.6, 0.8];

color = {soft_blue, warm_red, muted_green, soft_purple};

L2_v = zeros(2, length(Courant));
L2_p = zeros(2, length(Courant));

for i = 1:length(quad_type)
    for j = 1:length(Courant)
        disp(['using Cr: ', num2str(Courant(j))] );
        [L2_v(i, j) , L2_p(i, j)] = acoustic_solver_integrate(Courant(j), quad_type(i)); 
    end
end

font_size = 18;
figure(1)
semilogy(Courant, L2_p(1, :), '-', 'Color', warm_red ,'LineWidth', 1.5, 'MarkerFaceColor', warm_red, 'DisplayName', 'Gauss')
hold on
semilogy(Courant, L2_p(2, :), '-', 'Color', soft_blue ,'LineWidth', 1, 'MarkerFaceColor', soft_blue, 'DisplayName', 'Gauss Lobatto')
xlabel('Courant number (Cr)', 'FontSize', font_size)
ylabel('L_2 Pressure error', 'FontSize', font_size)
title('Change of L_2 Pressure Error vs. Courant Number', 'FontSize', font_size+2);
xlim([1, max(Courant)]);  
ylim([1e-13, 1e+1]);
set(gca, 'FontSize', font_size)
legend('FontSize', font_size)
hold off;

figure(2)
semilogy(Courant, L2_p(1, :), '-', 'Color', warm_red ,'LineWidth', 1.5, 'MarkerFaceColor', warm_red, 'DisplayName', 'Gauss')
hold on
semilogy(Courant, L2_p(2, :), '-', 'Color', soft_blue ,'LineWidth', 1, 'MarkerFaceColor', soft_blue, 'DisplayName', 'Gauss Lobatto')
xlabel('Courant number (Cr)', 'FontSize', font_size)
ylabel('L_2 Velocity error', 'FontSize', font_size)
title('Change of L_2 Velocity Error vs. Courant Number', 'FontSize', font_size+2);
xlim([1, max(Courant)]);  
ylim([1e-13, 1e+1]);
set(gca, 'FontSize', font_size)
legend('FontSize', font_size)
hold off;


function [l2_error_v, l2_error_p] = acoustic_solver_integrate(Courant, quad_type)
% 1D acoustic solver based on pointwise evaluation of fluxes and using
% integrals for the element advection term
% Assumption: Nodal polynomials with node points at interval end points
clear
close all

n = 80;         % number of elements
Tf = 2;        % final time
bc = 0;         % select boundary condition: Dirichlet (0), Absorbing(1)
k = 4;          % polynomial degree
c = 1;          % advection speed
rho = 1;        % density
Cr = Courant/k^2;   % Courant number -> sets time step size in dt = Cr * h / a
alpha = 0.0;    % flux type, 0 = upwind, 1 = central
nc = k+1;       % number of quadrature points
left = 0;       % left end of the domain
right = 1;      % right end of the domain
plot_accurate = 1; % plot only on nodes (0) or with more resolution (1)
flux_type = 0;  % choose between Lax-Friedrich (0) or HDG (1)

% analytical solution
analytical_v = @(x,t)cos(pi*x)*cos(pi*t);
analytical_p = @(x,t)sin(pi*x)*sin(pi*t);

% set quadrature formula and quadrature nodes for integration
if quad_type == 0
    [pg,wg] = get_gauss_quadrature(nc);
else
    [pg,wg] = get_gauss_lobatto_quadrature(nc);
end
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

% plot velocity numerical solution (red), the analytical solution (blue), and
% the initial condition
% figure(1)
% if plot_accurate == 1
%     xx_unit = -1:0.05:1;
%     xx = zeros(n,length(xx_unit));
%     vv_0 = zeros(size(xx));
%     vv = zeros(size(xx));
%     val = evaluate_lagrange_basis(xunit, xx_unit);
%     for e=1:n
%         xx(e, :) = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*xx_unit);
%         vv_0(e, :) = val' * w(e*kp1-k:e*kp1, 1);
%         vv(e, :) = val' * w(e*kp1-k:e*kp1, end);
%     end
%     v_anal = analytical_v(xx, Tf);
%     plot(xx(1,:),vv_0(1,:),'k:',xx(1,:),vv(1,:),'r-',xx(1,:),v_anal(1,:),'b');
%     hold on
%     plot(xx',vv_0','k:');
%     plot(xx',vv','r*-');
%     plot(xx',v_anal','b-');
%     hold off
% else
%     v_s = w(e*kp1-k:e*kp1, 1);
%     v_f = w(e*kp1-k:e*kp1, end);
%     plot(x,v_s,'k:', x, v_f,'r-', x, analytical_v(x, Tf),'b')
% end
% xlabel('x')
% ylabel('v_h(x)')
% title(['degree=' num2str(k) ', n=' num2str(n) ' elements, dt = ' num2str(dt)])
% legend('v_h(x,0)',['v_h(x,' num2str(Tf) ')'],['v(x,' num2str(Tf) ')'])
% 
% % plot pressure numerical solution (red), the analytical solution (blue), and
% % the initial condition
% figure(2)
% if plot_accurate == 1
%     xx_unit = -1:0.05:1;
%     xx = zeros(n,length(xx_unit));
%     pp_0 = zeros(size(xx));
%     pp = zeros(size(xx));
%     val = evaluate_lagrange_basis(xunit, xx_unit);
%     for e=1:n
%         xx(e, :) = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*xx_unit);
%         pp_0(e, :) = val' * w((e+n)*kp1-k: (e+n)*kp1, 1);
%         pp(e, :) = val' * w((e+n)*kp1-k: (e+n)*kp1, end);
%     end  
%     p_anal = analytical_p(xx, Tf);
%     plot(xx(1,:),pp_0(1,:),'k:',xx(1,:),pp(1,:),'r-',xx(1,:),p_anal(1,:),'b');
%     hold on
%     plot(xx',pp_0','k:');
%     plot(xx',pp','r-');
%     plot(xx',p_anal','b-');
%     hold off
% else
%     p_s = w((e+n)*kp1-k: (e+n)*kp1, 1);
%     p_f = w((e+n)*kp1-l: (e+n)*kp1, end);
%     plot(x,p_s,'k:', x, p_f,'r-', x, analytical_p(x, Tf),'b')
% end
% xlabel('x')
% ylabel('p_h(x)')
% title(['degree=' num2str(k) ', n=' num2str(n) ' elements, dt = ' num2str(dt)])
% legend('p_h(x,0)',['p_h(x,' num2str(Tf) ')'],['p(x,' num2str(Tf) ')'])

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
l2_error_v = l2error_v;

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
l2_error_p = l2error_p;

disp(['Pressure Error in maximum norm ' num2str(linfty_error_p) ' in L2 norm ' num2str(l2error_p)])

toc;
end

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

    % interpolate v and p to quadrature points
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
            case 1 % Absorbing
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
            case 1 % Absorbing
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
            case 1 % Absorbing check n
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
            case 1 % Absorbing check n
                lambda = (1 /(tau + 1/c)) * (rho * vminus * norm + tau * pminus);
        end
    else
        vplus = w(e*kp1+1); % e < n: w(kp1:(e-1)*kp1)
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