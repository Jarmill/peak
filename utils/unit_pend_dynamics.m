function [xdot] = unit_pend_dynamics(t, x, b, u)
%UNIT_PEND_DYNAMICS dynamics of a pendulum with all unit constants
%
%Input:
%   t:     time
%   x(1):  theta
%   x(2):  theta' = omega
%   b:     friction
%   u:     control law x -> u(x)

if nargin < 3
    b = 0;
end

if nargin < 4
    u = @(x) 0;
end

%the cosine is for cart-pole (?)
% xdot = [x(2); -sin(x(1)) + cos(x(1)) * u(x) - b*x(2)];
xdot = [x(2); -sin(x(1)) +  u(x) - b*x(2)];

end

