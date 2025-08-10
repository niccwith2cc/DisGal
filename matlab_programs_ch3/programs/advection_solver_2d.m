function advection_solver_2d
% 2D advection solver based on pointwise evaluation of fluxes and using
% integrals for the element advection term, using a variable transport
% speed
% Assumption: Nodal polynomials with node points at interval end points

N = [20, 20];    % number of elements
Tf = 2;         % final time
k = 5;          % polynomial degree
Cr = 0.5/k^1.5; % Courant number -> sets time step size in dt = Cr * h / a
alpha = 0.0;    % flux type, 0 = upwind, 1 = central
nc = k+1;       % number of quadrature points
domain = [0 1 0 1]; % domain as rectangle [xmin xmax ymin ymax]
output_freq = 0.1;  % how often to print an output figure

% variable advection speed
b = @(x,y,t)cos(pi/Tf*t)*[2*sin(2*pi*y).*sin(pi*x).^2,...
                          -2*sin(2*pi*x).*sin(pi*y).^2];
%b = @(x,y,t)[ones(size(x)) ones(size(x))];

analytical = @(x,y,t)exp(-400*((x-0.5).^2+(y-0.75).^2));
%analytical = @(x,y,t)sin(2*pi*(y-t));

% set 1D quadrature formula and quadrature nodes for integration
[pg,wg] = get_gauss_quadrature(nc); 

% set the node points for the Lagrange polynomials
xunit = get_gauss_lobatto_quadrature(k+1);

% start up the simulation

% create mesh y compute mesh size h
n = N(1) * N(2);
h = [(domain(2)-domain(1))/N(1), (domain(4)-domain(3))/N(2)];

% compute time step from Cr number with max velocity 2, adjust to hit the
% desired final time exactly
dt = Cr * min(h) / 2;
NT = round(Tf/dt);
dt = Tf/NT;

disp(['Number of elements: ' num2str(n) ', number of DoFs: ' num2str(n*(k+1)^2)])
disp([' minimum mesh size: ' num2str(min(h)) ', time step size: ' num2str(dt) ])

% evaluate reference cell polynomials and mass matrix
[values,derivatives] = evaluate_lagrange_basis(xunit, pg);
Me = 0.25 * h(1) * h(2) * ...
    (kron(values,values) * diag(kron(wg,wg)) * kron(values,values)');
invMe = inv(Me);
Minv = sparse((k+1)^2*n,(k+1)^2*n);
for e=1:n
    Minv((e-1)*(k+1)^2+1:e*(k+1)^2,(e-1)*(k+1)^2+1:e*(k+1)^2) = invMe;
end

tic

% set initial condition by projection
u=zeros((k+1)^2*n,1);
for e_y=1:N(2)
    for e_x=1:N(1)
        e = e_x + N(1)*(e_y-1);
        % compute points in reference domain
        x = kron(ones(nc,1), domain(1)+(e_x-1+0.5+0.5*pg)*h(1));
        y = kron(domain(3)+(e_y-1+0.5+0.5*pg)*h(2), ones(nc,1));
        u((e-1)*(k+1)^2+1:e*(k+1)^2) = kron(values,values)'\analytical(x,y,0);
    end
end

[mass,energy] = compute_balance(domain, k, N, xunit, u);
disp(['t=0.000 mass ' num2str(mass,12) ' energy ' num2str(energy,12)])

figure(1)
u_init = u;
plot_solution(domain, k, N, xunit, 0, u)

% run time loop
for m=1:NT
    k1 = Minv*evaluate_adve_rhs(u, b, domain, N, (m-1)*dt, ...
                 analytical, values, derivatives, pg, wg, alpha);
    k2 = Minv*evaluate_adve_rhs(u+0.5*dt*k1, b, domain, N, (m-0.5)*dt, ...
                 analytical, values, derivatives, pg, wg, alpha);
    k3 = Minv*evaluate_adve_rhs(u+0.5*dt*k2, b, domain, N, (m-0.5)*dt, ...
                 analytical, values, derivatives, pg, wg, alpha);
    k4 = Minv*evaluate_adve_rhs(u+dt*k3, b, domain, N, m*dt, ...
                 analytical, values, derivatives, pg, wg, alpha);
    
    u = u + dt/6*(k1+2*k2+2*k3+k4);
    
    if floor(m*dt / output_freq) ~= floor((m-1)*dt / output_freq)
        [mass,energy] = compute_balance(domain, k, N, xunit, u);
        disp(['t=' num2str(m*dt,'%5.3f') ' mass ' num2str(mass,12) ...
            ' energy ' num2str(energy,12)])
        plot_solution(domain, k, N, xunit, m*dt, u)
    end
end

toc

disp(['Distance of node values ' num2str(max(abs(u-u_init)))])

end


function rhs = evaluate_adve_rhs(u, b, domain, N, t, bc, values, derivatives, pg, weights, alpha)
% evaluate the right hand side of the two-dimensional advection equation on
% a rectangular grid, defined by the domain, number of points in the two
% directions and one-dimensional interpolation matrices. The
% multi-dimensional interpolation matrices are defined by the tensor
% product (Kronecker product) of the 1D matrices

kp1 = size(values, 1); % degree + 1
k = kp1 - 1;
nc = size(values, 2);
h = [(domain(2)-domain(1))/N(1)  (domain(4)-domain(3))/N(2)];
rhs = zeros(size(u));

% two-dimensional interpolation matrices and quadrature weights by
% Kronecker product
values_2d = kron(values,values);
derivatives_2d = [kron(values, derivatives) kron(derivatives, values)];
weights_2d = kron(weights,weights);

% indices of the nodes at the faces
face_indices = zeros(kp1, 4);
face_indices(:,1) = 1:kp1:k*kp1+1; % left
face_indices(:,2) = kp1:kp1:kp1^2; % right
face_indices(:,3) = 1:kp1;         % bottom
face_indices(:,4) = k*kp1+1:kp1^2; % top
neighbor_face_indices=[face_indices(:,2),face_indices(:,1),...
                       face_indices(:,4),face_indices(:,3)];
xbound = [-ones(nc,1), ones(nc,1), pg, pg];
ybound = [pg, pg, -ones(nc,1) ones(nc,1)];
ref_x = kron(ones(nc,1), pg);
ref_y = kron(pg, ones(nc,1));
% tangentials oriented in clockwise direction
ref_tangentials = [0 0 -1 1; 1 -1 0 0];
neighbor_stride = [-1 0; 1 0 ; 0 -1; 0 1];

% we have the same Jacobian and determinant on all elements, so do this
% outside of the loop
jac = [(domain(2)-domain(1))/N(1)/2 0 ; 0 (domain(4)-domain(3))/N(2)/2];
jacdet = det(jac);
invjac = inv(jac);

for e_y=1:N(2)
    for e_x=1:N(1)
        e = e_x + N(1)*(e_y-1);
        
        % compute velocity at quadrature points
        x = domain(1) + (e_x-0.5+0.5*ref_x)*h(1);
        y = domain(3) + (e_y-0.5+0.5*ref_y)*h(2);
        be = b(x,y,t);
        
        % extract element solution
        ue = u((e-1)*kp1^2+1:e*kp1^2);
        
        % interpolate u to quadrature points
        u_quad = (values_2d')*ue;
               
        % multiply be by inverse Jacobian
        j_be = be * invjac';
        % the line above here is equivalent to
        % j_be = [be(:,1)*invjac(1,1)+be(:,2)*invjac(1,2) ...
        %         be(:,1)*invjac(2,1)+be(:,2)*invjac(2,2)];
                
        % compute operator at quadrature points and multiply by gradient
        % of test function
        x_f = j_be .* weights_2d .* jacdet .* u_quad;
        myrhs = derivatives_2d * x_f(:);

        % comment in this evaluation of the function on all faces at once 
        % to get better performance
        %bf_all = b(domain(1)+(e_x-0.5+0.5*xbound)*h(1),...
        %           domain(3)+(e_y-0.5+0.5*ybound)*h(2),t);
    
        % faces
        for f=1:4
            % interpolate interior solution to quadrature points
            uminus = values'*ue(face_indices(:,f));

            % resolve boundary values or interpolate neighboring solution
            % to quadrature points
            if (f == 1 && e_x == 1) || (f == 2 && e_x == N(1)) || ...
               (f == 3 && e_y == 1) || (f == 4 && e_y == N(2))
                uplus = bc(domain(1)+(e_x-0.5+0.5*xbound(:,f))*h(1),...
                           domain(3)+(e_y-0.5+0.5*ybound(:,f))*h(2),t);
            else
                eplus = neighbor_stride(f,:) * [1; N(1)] + e;
                ueplus = u((eplus-1)*kp1^2+1:eplus*kp1^2);
                uplus = values'*ueplus(neighbor_face_indices(:,f));
            end
            tangent = jac * ref_tangentials(:,f);
            scale = norm(tangent);
            normal = [-tangent(2); tangent(1)]/scale;
            
            % faster: comment in code above and this line here and comment
            % out the next line
            % bf = [bf_all(:,f) bf_all(:,4+f)];
            bf = b(domain(1)+(e_x-0.5+0.5*xbound(:,f))*h(1),...
                   domain(3)+(e_y-0.5+0.5*ybound(:,f))*h(2),t);
            
            % compute numerical flux at all quadrature points at once and
            % multiply it by the normal vector
            numflux_n = (uminus.*bf+uplus.*bf)/2*normal + ...
                (1-alpha)/2*abs(bf * normal) .* (uminus - uplus);
                
            % multiply by test functions, sum over quadrature points
            myrhs(face_indices(:,f)) = myrhs(face_indices(:,f)) - ...
                values * (scale * weights .* numflux_n);
        end
        
        rhs((e-1)*kp1^2+1:e*kp1^2) = myrhs;
    end
end

end

function plot_solution(domain, k, N, xunit, t, u)
% function to plot the solution represented in an element-wise way with the
% surf function, i.e., resort the unknowns in a matrix with lexicographic
% order and call surf()

h = [(domain(2)-domain(1))/N(1)  (domain(4)-domain(3))/N(2)];
n = N(1) * N(2);
% prepare indices
xplot = zeros(N(1)*(k+1),1);
for e=1:N(1)
    xplot(e*(k+1)-k:e*(k+1)) = domain(1)+(e-1+0.5+0.5*xunit)*h(1);
end
yplot = zeros(N(2)*(k+1),1);
for e=1:N(2)
    yplot(e*(k+1)-k:e*(k+1)) = domain(3)+(e-1+0.5+0.5*xunit)*h(2);
end
ele2lex = zeros(n*(k+1)^2,1);
for e_y=1:N(2)
    for e_x=1:N(1)
        for j=0:k
            start = 1+((e_y-1)*(k+1)+j)*N(1)*(k+1) + (e_x-1)*(k+1);
            starti = 1+((e_y-1)*N(1)+(e_x-1))*(k+1)^2 + j*(k+1);
            ele2lex(start:start+k) = starti:starti+k;
        end
    end
end

surf(xplot,yplot,reshape(u(ele2lex),N(1)*(k+1),N(2)*(k+1))');
xlabel('x')
ylabel('y')
zlabel(['u(x,y,' num2str(t) ')'])
axis([domain -0.2 1])
pause(0.05)

end

function [mass, energy] = compute_balance(domain, k, N, xunit, u)

[pg_err,wg_err] = get_gauss_quadrature(k+1);
values_err = evaluate_lagrange_basis(xunit, pg_err);
values_err_2d = kron(values_err, values_err);
wg_err_2d = kron(wg_err, wg_err);

jacdet = (domain(2)-domain(1))/N(1) * (domain(4)-domain(3))/N(2) * 0.25;
n = N(1) * N(2);

mass = 0;
energy = 0;
for e=1:n
    sol_num = values_err_2d' * u((e-1)*(k+1)^2+1:e*(k+1)^2);
    mass = mass + jacdet * wg_err_2d' * sol_num;
    energy = energy + jacdet * wg_err_2d' * sol_num.^2;
end

end
