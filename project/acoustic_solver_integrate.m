% function acoustic_solver_integrate
% 1D acoustic solver based on pointwise evaluation of fluxes and using
% integrals for the element advection term
% Assumption: Nodal polynomials with node points at interval end points
close all
clear

n = 20;         % number of elements
Tf = 2;         % final time
periodic = 0;   % switch between Dirichlet conditions (0) and periodic (1)
k = 2;          % polynomial degree
c = 1;         % advection speed
rho = 1;       % density
Cr = 0.4/k^2;   % Courant number -> sets time step size in dt = Cr * h / a
alpha = 0.0;    % flux type, 0 = upwind, 1 = central
nc = k+1;       % number of quadrature points
left = 0;       % left end of the domain
right = 1;      % right end of the domain
plot_accurate = 1; % plot only on nodes (0) or with more resolution (1)

% analytical solution
% analytical = @(x,t)sin(4*pi*(x-a*t));
% analytical = @(x,t)exp(sin(4*pi*(x-a*t)));
% analytical = @(x,t)(abs(2*mod(x-a*t,1)-1));
% analytical = @(x,t)(mod(x-a*t,1)>0.5);
analytical_v = @(x,t)cos(pi*x)*cos(pi*t);
analytical_p = @(x,t)sin(pi*x)*sin(pi*t);


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

disp(['Number of elements: ' num2str(n) ', minimum mesh size: ' ...
    num2str(min(h)) ', time step size: ' num2str(dt) ])

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
    k1 = Minv * evaluate_acoustic_rhs(w(:, m), c, rho, [analytical_v(0, (m-1)*dt); analytical_v(1, (m-1)*dt); analytical_p(0, (m-1)*dt); analytical_p(1, (m-1)*dt)], values, derivatives, wg, alpha, periodic);

    k2 = Minv * evaluate_acoustic_rhs(w(:, m) + 0.5*dt*k1, c, rho, [analytical_v(0, (m-0.5)*dt); analytical_v(1, (m-0.5)*dt); analytical_p(0, (m-0.5)*dt); analytical_p(1, (m-0.5)*dt)], values, derivatives, wg, alpha, periodic);
                 
    k3 = Minv * evaluate_acoustic_rhs(w(:, m) + 0.5*dt*k2, c, rho, [analytical_v(0, (m-0.5)*dt) ; analytical_v(1, (m-0.5)*dt); analytical_p(0, (m-0.5)*dt) ; analytical_p(1, (m-0.5)*dt)], values, derivatives, wg, alpha, periodic);

    k4 = Minv * evaluate_acoustic_rhs(w(:, m) + dt*k3, c, rho, [analytical_v(0, m*dt); analytical_v(1, m*dt); analytical_p(0, m*dt); analytical_p(1, m*dt)], values, derivatives, wg, alpha, periodic);
    
    w(:, m+1) = w(:, m) + dt/6*(k1+2*k2+2*k3+k4);
end

% plot the numerical solution (red), the analytical solution (blue), and
% the initial condition
figure(1)
if plot_accurate == 1
    xx_unit = -1:0.05:1;
    xx = zeros(n,length(xx_unit));
    vv_0 = zeros(size(xx));
    vv = zeros(size(xx));
    val = evaluate_lagrange_basis(xunit, xx_unit);
    for e=1:n
        xx(e, :) = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*xx_unit);
        vv_0(e, :) = val' * w(e*kp1-k:e*kp1, 1);
        vv(e, :) = val' * w(e*kp1-k:e*kp1, end);
    end
    v_anal = analytical_v(xx, Tf);
    plot(xx(1,:),vv_0(1,:),'k:',xx(1,:),vv(1,:),'r-',xx(1,:),v_anal(1,:),'b');
    hold on
    plot(xx',vv_0','k:');
    plot(xx',vv','r-');
    plot(xx',v_anal','b-');
    hold off
else
    v_s = w(e*kp1-k:e*kp1, 1);
    v_f = w(e*kp1-k:e*kp1, end);
    plot(x,v_s,'k:', x, v_f,'r-', x, analytical_v(x, Tf),'b')
end
xlabel('x')
ylabel('v_h(x)')
title(['degree=' num2str(k) ', n=' num2str(n) ' elements, dt = ' num2str(dt)])
legend('v_h(x,0)',['v_h(x,' num2str(Tf) ')'],['v(x,' num2str(Tf) ')'])

figure(2)
if plot_accurate == 1
    xx_unit = -1:0.05:1;
    xx = zeros(n,length(xx_unit));
    pp_0 = zeros(size(xx));
    pp = zeros(size(xx));
    val = evaluate_lagrange_basis(xunit, xx_unit);
    for e=1:n
        xx(e, :) = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*xx_unit);
        pp_0(e, :) = val' * w((e+n)*kp1-k: (e+n)*kp1, 1);
        pp(e, :) = val' * w((e+n)*kp1-k: (e+n)*kp1, end);
    end  
    p_anal = analytical_p(xx, Tf);
    plot(xx(1,:),pp_0(1,:),'k:',xx(1,:),pp(1,:),'r-',xx(1,:),p_anal(1,:),'b');
    hold on
    plot(xx',pp_0','k:');
    plot(xx',pp','r-');
    plot(xx',p_anal','b-');
    hold off
else
    p_s = w((e+n)*kp1-k: (e+n)*kp1, 1);
    p_f = w((e+n)*kp1-l: (e+n)*kp1, end);
    plot(x,p_s,'k:', x, p_f,'r-', x, analytical_p(x, Tf),'b')
end
xlabel('x')
ylabel('p_h(x)')
title(['degree=' num2str(k) ', n=' num2str(n) ' elements, dt = ' num2str(dt)])
legend('p_h(x,0)',['p_h(x,' num2str(Tf) ')'],['p(x,' num2str(Tf) ')'])

l2error = 0;
linfty_error = 0;
[pg_err,wg_err] = get_gauss_quadrature(k+3);
values_err = evaluate_lagrange_basis(xunit, pg_err);
for e=1:n
    sol_num = values_err' * w((e+n-1)*kp1+1:(e+n)*kp1, end);
    x_err = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*pg_err);
    sol_exact = analytical_p(x_err, Tf);
    l2error = l2error + h(e)/2 * wg_err' * (sol_num-sol_exact).^2;
    linfty_error = max([linfty_error; abs(sol_num-sol_exact)]);
end
l2error = sqrt(l2error);

disp(['Error in maximum norm ' num2str(linfty_error) ' in L2 norm ' num2str(l2error)])
toc;

% end


function rhs = evaluate_acoustic_rhs(w, c, rho, bc, values, derivatives, weights, alpha, periodic)
% disp(['this is velocity bc: ', num2str(bc(1))])

kp1 = size(values, 1); % degree + 1
n = length(w)/(2*kp1);
rhs = zeros(size(w));

for e=1:n
    ve = w((e-1)*kp1+1: e*kp1); % 1 to kp1*n
    pe = w((e+n-1)*kp1+1: (e+n)*kp1); % from kp1*n+1 to 2*kp1*n
    disp(['size of pe: ', num2str(length(pe))])

    % interpolate v to quadrature points
    v_quad = values' * ve;
    p_quad = values' * pe;

    flux_v = weights .* p_quad;
    flux_p = weights .* v_quad;
    % compute operator at quadrature points and multiply by gradient of
    % test function
    rhs((e-1)*kp1+1: e*kp1) = derivatives * flux_v;
    rhs((e+n-1)*kp1+1: (e+n)*kp1) = derivatives * flux_p;
    
    % compute acoustic numerical flux on the left
    vminus = ve(1);
    pminus = pe(1);
    if (e==1)
        if (periodic)
            % periodic bc
            vplus = w(kp1*n); % last value in the velocity
            pplus = w(2*kp1*n);
        else
            % Dirichlet condition, implemented via mirror principle
%             vplus = 2*bc(1)-vminus;
            pplus = 2*bc(3)-pminus;
            vplus = -ve(1);
        end
    else
        vplus = w((e-1)*kp1); % e > 1: w(kp1:(e-1)*kp1)
        pplus = w((e+n-1)*kp1); 
    end
    numflux_v = 1/2*c*rho^2 * (vminus + vplus) + 1/2 * abs(c) * (pminus - pplus);
    numflux_p = 1/(2*rho) * (pminus + pplus) + 1/2 * abs(c) * (vminus - vplus);
    rhs((e-1)*kp1+1) = rhs((e-1)*kp1+1) + numflux_v;
    rhs((e+n-1)*kp1+1) = rhs((e+n-1)*kp1+1) + numflux_p;
    
    % compute acoustic numerical flux on the right
    vminus = ve(kp1);
    pminus = pe(kp1);
    if (e==n)
        if (periodic)
            % periodic bc
            vplus = w(1);
            pplus = w(kp1*n+1);
        else
            % Dirichlet condition, implemented via mirror principle
%             vplus = 2*bc(2)-vminus;
            pplus = 2*bc(4)-pminus;
            vplus = -ve(end);
        end
    else
        vplus = w(e*kp1+1);
        pplus = w((e+n)*kp1+1);
    end
    numflux_v = 1/2*c*rho^2 * (vminus + vplus) + 1/2 * abs(c)*(pplus - pminus);
    numflux_p = 1/(2*rho) * (pminus + pplus) + 1/2 * abs(c)*(vplus - vminus);
    rhs(e*kp1) = rhs(e*kp1) - numflux_v;
    rhs((e+n)*kp1) = rhs((e+n)*kp1) - numflux_p;
end

end

% function rhs = evaluate_acoustic_rhs_hdg(w, c, rho, bc, values, derivatives, weights, tau, periodic)
% % HDG flux implementation for 1D acoustic wave equation
% % Inputs:
% %   w        - full solution vector: [v; p]
% %   c, rho   - physical constants
% %   bc       - [vL, vR, pL, pR] Dirichlet boundary values
% %   tau      - stabilization parameter (e.g., 0.5 or 1.0)
% 
% kp1 = size(values, 1); % degree + 1
% n = length(w)/(2*kp1);
% rhs = zeros(size(w));
% 
% for e = 1:n
%     ve = w((e-1)*kp1+1: e*kp1); 
%     pe = w((e+n-1)*kp1+1: (e+n)*kp1); 
% 
%     % interpolate to quadrature points
%     v_quad = values' * ve;
%     p_quad = values' * pe;
% 
%     % physical terms
%     rhs((e-1)*kp1+1: e*kp1) = derivatives * (weights .* p_quad);
%     rhs((e+n-1)*kp1+1: (e+n)*kp1) = derivatives * (weights .* v_quad);
%     
%     % LEFT FACE
%     vminus = ve(1);
%     pminus = pe(1);
%     if e == 1
%         if periodic
%             vplus = w(kp1*n);
%             pplus = w(2*kp1*n);
%         else
%             pplus = bc(3);  % Dirichlet p
%             vplus = bc(1);  % Dirichlet v
%         end
%     else
%         vplus = w((e-1)*kp1);
%         pplus = w((e+n-1)*kp1);
%     end
% 
%     avg_p = 0.5 * (pminus + pplus);
%     jump_v = vplus - vminus;
%     flux_v = avg_p - tau * jump_v;
% 
%     avg_v = 0.5 * (vminus + vplus);
%     jump_p = pplus - pminus;
%     flux_p = avg_v - tau * jump_p;
% 
%     rhs((e-1)*kp1+1) = rhs((e-1)*kp1+1) + flux_v;
%     rhs((e+n-1)*kp1+1) = rhs((e+n-1)*kp1+1) + flux_p;
% 
%     % RIGHT FACE
%     vminus = ve(end);
%     pminus = pe(end);
%     if e == n
%         if periodic
%             vplus = w(1);
%             pplus = w(kp1*n+1);
%         else
%             pplus = bc(4);  % Dirichlet p
%             vplus = bc(2);  % Dirichlet v
%         end
%     else
%         vplus = w(e*kp1+1);
%         pplus = w((e+n)*kp1+1);
%     end
% 
%     avg_p = 0.5 * (pminus + pplus);
%     jump_v = vplus - vminus;
%     flux_v = avg_p - tau * jump_v;
% 
%     avg_v = 0.5 * (vminus + vplus);
%     jump_p = pplus - pminus;
%     flux_p = avg_v - tau * jump_p;
% 
%     rhs(e*kp1) = rhs(e*kp1) - flux_v;
%     rhs((e+n)*kp1) = rhs((e+n)*kp1) - flux_p;
% end
% end

