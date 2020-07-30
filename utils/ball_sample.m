function [X] = ball_sample(N, d)
%BALL_SAMPLE Randomly sample a point from a unit d-dimensional ball

%dropped coordinate method
%U = randn(N, d+2);
%normU = sqrt(sum(U.^2, 2));
%u = U ./ normU;
u = sphere_sample(N, d+2);
X = u(:, 1:d);

end

