function [points, weights] = get_gauss_lobatto_quadrature(n_points)
% get the position and the weights of the Gauss-Lobatto quadrature with the
% given number of points

points = zeros(n_points, 1);
weights = zeros(n_points, 1);
if n_points > 2
    M=zeros(n_points-2,n_points-2);
    for i=1:n_points-3
        M(i,i+1) = sqrt(4*i*(i+1)*(i+1)*(i+2)/((2*i+1)*(2*i+2)^2*(2*i+3)));
        M(i+1,i) = M(i,i+1);
    end
    eig(M);
    points(2:n_points-1) = eig(M);
end
points(1) = -1;
points(n_points) = 1;
for i=1:n_points
    weights(i) = 2*factorial(n_points-1)^2/((n_points-1)*factorial(n_points-1)*...
        factorial(n_points)*jacobi_polynomial(points(i), 0, 0, n_points-1)^2);
end
end


function val = jacobi_polynomial(x, alpha, beta, n)

p0 = 1;
p1 = ((alpha+beta+2)*x + (alpha-beta))/2;
if n==0
    val = p0;
elseif n==1
    val = p1;
else
    for i=1:n-1
        v = 2*i+(alpha+beta);
        a1 = 2*(i+1)*(i+alpha+beta+1)*v;
        a2 = (v+1)*(alpha^2-beta^2);
        a3 = v*(v+1)*(v+2);
        a4 = 2*(i+alpha)*(i+beta)*(v+2);
        pn = ((a2+a3*x)*p1-a4*p0)/a1;
        p0 = p1;
        p1 = pn;
    end
    val = p1;
end
end