function [X] = circle_sample(N)
%CIRCLE_SAMPLE Uniformly Sample N points from the unit circle

    theta = 2*pi*rand(N, 1);
    trig = [cos(theta) sin(theta)];
    r = sqrt(rand(N, 1));

    X = r.*trig;
end

