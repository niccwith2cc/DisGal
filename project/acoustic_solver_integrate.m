function acoustic_solver_integrate
% 1D acoustic solver based on pointwise evaluation of fluxes and using
% integrals for the element advection term
% Assumption: Nodal polynomials with node points at interval end points
close all
clear

n = 40;         % number of elements
Tf = 2;         % final time
periodic = 0;   % switch between Dirichlet conditions (0) and periodic (1)
k = 5;          % polynomial degree
c = +1;         % advection speed
rho = +1;       % density
Cr = 0.4/k^2;   % Courant number -> sets time step size in dt = Cr * h / a
alpha = 0.0;    % flux type, 0 = upwind, 1 = central
nc = k+1;       % number of quadrature points
left = 0;       % left end of the domain
right = 1;      % right end of the domain
plot_accurate = 1; % plot only on nodes (0) or with more resolution (1)

% analytical solution
%analytical = @(x,t)sin(4*pi*(x-a*t));
% analytical = @(x,t)exp(sin(4*pi*(x-a*t)));
% analytical = @(x,t)(abs(2*mod(x-a*t,1)-1));
% analytical = @(x,t)(mod(x-a*t,1)>0.5);
analytical_p = @(x,t)sin(pi*x).*sin(pi*t);
analytical_v = @(x,t)cos(pi*x).*cos(pi*t);

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
Minv = sparse(2*kp1*n,2*kp1*n);
for e=1:n
    Mloc = blkdiag(inv(0.5*h(e)*Me), inv(0.5*h(e)*Me));
    Minv((2*kp1)*(e-1)+ 1:(2*kp1)*e, (2*kp1)*(e-1)+1:(2*kp1)*e) = Mloc;
end

% set initial condition
v0=@(x) cos(pi*x);
p0=@(x) sin(pi*x);
w0 = [v0(x); p0(x)];
w = zeros(2*kp1*n, NT+1);
w(:, 1) = w0;

% run time loop
for m=1:NT
    k1 = Minv * evaluate_acoustic_rhs(w(:,m), values, derivatives, wg, alpha, periodic);

    k2 = Minv * evaluate_acoustic_rhs(w(:,m) + 0.5*dt*k1, values, derivatives, wg, alpha, periodic);
                 
    k3 = Minv * evaluate_acoustic_rhs(w(:,m) + 0.5*dt*k2, values, derivatives, wg, alpha, periodic);

    k4 = Minv * evaluate_acoustic_rhs(w(:,m) + dt*k3, values, derivatives, wg, alpha, periodic);
    
    w(:,m+1) = w(:,m) + dt/6*(k1+2*k2+2*k3+k4);
end

% plot the numerical solution (red), the analytical solution (blue), and
% the initial condition
figure(1)
if plot_accurate == 1
    for e=1:n
        xx_unit = linspace(y(2*e-1), y(2*e), 50);
        [~, val] = evaluate_lagrange_basis(xunit, 2*(xx_unit - y(2*e-1))/h(e) - 1);
        we = w((2*kp1)*(e-1)+1:(2*kp1)*e, end);
        v_num = val' * we(1:kp1);
        p_num = val' * we(kp1+1:end);
        plot(xx_unit, p_num, 'r-', xx_unit, p_anal(xx_unit, Tf), 'b--'); hold on;
    end
    legend("Numerical p", "Exact p");
    xlabel("x"); ylabel("(p(x, Tf)");
    title("1D Acoustic Wave Equation");
else
%     plot(x,u(:,1),'k:',x,u(:,end),'r-',x,analytical(x,Tf),'b')
end
% xlabel('x')
% ylabel('u_h(x)')
% title(['degree=' num2str(k) ', n=' num2str(n) ' elements, dt = ' num2str(dt)])
% legend('u_h(x,0)',['u_h(x,' num2str(Tf) ')'],['u(x,' num2str(Tf) ')'])
% 
% l2error = 0;
% linfty_error = 0;
% [pg_err,wg_err] = get_gauss_quadrature(k+3);
% values_err = evaluate_lagrange_basis(xunit, pg_err);
% for e=1:n
%     sol_num = values_err' * u((e-1)*kp1+1:e*kp1,end);
%     x_err = y(2*e-1)+(y(2*e)-y(2*e-1))*(0.5+0.5*pg_err);
%     sol_exact = analytical(x_err, Tf);
%     l2error = l2error + h(e)/2 * wg_err' * (sol_num-sol_exact).^2;
%     linfty_error = max([linfty_error; abs(sol_num-sol_exact)]);
% end
% l2error = sqrt(l2error);
%     
% disp(['Error in maximum norm ' num2str(linfty_error) ' in L2 norm ' num2str(l2error)])
toc;

end


function rhs = evaluate_acoustic_rhs(w, values, derivatives, wg, alpha, periodic)

kp1 = size(values, 1); % degree + 1
n = length(w)/kp1;
rhs = zeros(size(w));
c = 1; rho = 1;

for e = 1:n
    ind = (2*kp1)*(e-1)+1 : (2*kp1)*e;
    ve = w(ind(1:kp1));
    pe = w(ind(kp1+1:end));
    
    vq = values'*ve; pq = values'*pe;

    rhs_v = derivatives * ((diag(wg)) * pq);
    rhs_p = derivatives * (diag(wg) * vq);

    rhs(ind(1:kp1)) = rhs_v;
    rhs(ind(kp1+1:end)) = rhs_p;
    
    % Numerical flux
%     vminus = ve(1);
%     pminus = pe(1);
    if (e==1)
        if (periodic)
            % periodic bc
            vplus = w((2*kp1)*(n-1)+kp1);
            pplus = w((2*kp1)*(n-1)+ 2*kp1);
        else
            % Dirichlet condition, implemented via mirror principle
            vplus = ve(1);
            pplus = -pe(1);
        end
    else
        vplus = w((2*kp1)*(e-2)+kp1); 
        pplus = w((2*kp1)*(e-2)+2*kp1);
    end
    vminus = ve(1);
    pminus = pe(1);
    
    v_star = 0.5*(vminus + vplus) + (1-alpha)/(2*rho*c)*(pminus - pplus);
    p_star = 0.5*(pminus + pplus) + (1-alpha)*0.5*c*rho(vminus - vplus);

    rhs(ind(1)) = rhs(ind(1)) + (p_star - pminus);
    rhs(ind(kp1+1)) = rhs(ind(kp1+1)) + (v_star  - vminus);

    % compute advective numerical flux on the right
%     vminus = ve(kp1);
%     pminus = pe(kp1);
    if (e==n)
        if (periodic)
            % periodic bc
            vplus = w(kp1);
            pplus = w(2*kp1);
        else
            % Dirichlet condition, implemented via mirror principle
           vplus = ve(end);
           pplus = -pe(end);
        end
    else
        vplus = w((2*kp1)*e+1);
        pplus = w((2*kp1)*e+kp1+1);
    end
    vminus = ve(end);
    pminus = pe(end);

    v_star = 0.5*(vminus + vplus) + (1-alpha)/(2*rho*c)*(pminus - pplus);
    p_star = 0.5*(pminus + pplus) + (1-alpha)*0.5*c*rho(vminus - vplus);

    rhs(ind(kp1)) = rhs(ind(kp1)) - (p_star - pminus);
    rhs(ind(end)) = rhs(ind(end)) - (v_star  - vminus);  
end

end
