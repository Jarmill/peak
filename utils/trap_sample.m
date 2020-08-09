function [X] = trap_sample(N, I_max)
% sample from trapezoid [0<=x<=1], [0<=y<=I_max]

x = rand(15*N, 1);
y = rand(15*N, 1) * (I_max);

rej = (x+y) > 1;

x(rej) = [];
y(rej) = [];

X = [x(1:N), y(1:N)]';


end