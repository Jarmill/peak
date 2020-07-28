function [X] = sphere_sample(N, d)
%BALL_SAMPLE Randomly sample a point from a unit d-dimensional sphere
%3d sphere is a 2-sphere. S^(d-1)

%dropped coordinate method
U = randn(N, d);
normU = sqrt(sum(U.^2, 2));
X = U ./ normU;



end

