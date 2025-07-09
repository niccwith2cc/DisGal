function [values, derivatives] = evaluate_lagrange_basis(points_lagrange,points_evaluation)
% construct a Lagrange basis in the given points_lagrange and evaluate the
% Lagrange polynomials and their first derivative in a second of points of
% the argument points_evaluation

values = zeros(length(points_lagrange), length(points_evaluation));
derivatives = zeros(length(points_lagrange), length(points_evaluation));

% compute the Lagrange weights for all polynomials
weights = zeros(size(points_lagrange));
for i=1:length(weights)
    weights(i) = 1;
    for j=1:length(weights)
        if (j~=i)
            weights(i) = weights(i) * (points_lagrange(i)-points_lagrange(j));
        end
    end
    weights(i) = 1./weights(i);
end

% evaluation with the barycentric formula, see e.g. Berrut & Trefethen,
% SIAM Review 46(3):501-517, 2004
for q=1:length(points_evaluation)
    p = points_evaluation(q);
    for i=1:length(points_lagrange)
        value = weights(i);
        derivative = 0;
        for k=1:length(points_lagrange)
            if (i~=k)
                derivative = derivative * (p - points_lagrange(k)) + value;
                value = value * (p - points_lagrange(k));
            end
        end
        values(i,q) = value;
        derivatives(i,q) = derivative;
    end
end