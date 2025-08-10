function [points, weights] = get_gauss_quadrature(n_points)
% get the position and the weights of the Gaussian quadrature with the
% given number of points

points = zeros(n_points, 1);
weights = zeros(n_points, 1);
m = floor((n_points+1)/2);
for i=1:m
    z = cos(pi*(i-0.25)/(n_points+0.5));
    pp = 1;
    p1 = 1;
    % Newton iteration to find zeros of Jacobi polynomials with
    % alpha=beta=0 (Legendre polynomial)
    while (abs(p1/pp)>1e-13)
        p0 = 0;
        p1 = 1;
        % evaluate Legendre polynomial by recursion formula
        for j=1:n_points
            pold = p0;
            p0 = p1;
            p1 = ((2*j-1)*z*p0-(j-1)*pold)/j;
        end
        % derivative through Bonett's formula
        pp = n_points*(z*p1-p0)/(z^2-1);
        z = z-p1/pp;
    end
    points(i) = -z;
    points(n_points-i+1) = z;
    weights(i) = 2/((1-z^2)*pp^2);
    weights(n_points-i+1) = weights(i);
end


