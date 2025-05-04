function [values, derivatives] = evaluate_legendre_basis(degree,points_evaluation)
% construct a Lagrange basis in the given points_lagrange and evaluate the
% Lagrange polynomials and their first derivative in a second of points of
% the argument points_evaluation

values = zeros(degree+1, length(points_evaluation));
derivatives = zeros(degree+1, length(points_evaluation));

% Make sure evaluation works for both row and column vectors
x = ones(1,length(points_evaluation));
x(:) = points_evaluation;

values(1,:) = 1/sqrt(2);
if degree > 0
    values(2,:) = sqrt(1.5) * x;
    derivatives(2,:) = sqrt(1.5);
end
aj = sqrt(1/3);
% This is Bonnet's recursion formula with coefficients adapted to the
% scaling with unit matrix for the mass matrix
for j=2:degree
    ajm = aj;
    aj = sqrt(j^2/(2*j+1)/(2*j-1));
    values(j+1,:) = 1/aj * x .* values(j,:) - ajm / aj * values(j-1,:);
    derivatives(j+1,:) = 1/aj * (values(j,:) + x .* derivatives(j,:)) -...
        ajm / aj * derivatives(j-1,:);
end

end