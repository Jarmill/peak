% Test plotting a vector field
% Example from 'On Analysis and Synthesis of Safe Control Laws'
% by Anders Rantzer and Stephen Prajna

m = 4;
N = 20;

[x, y] = meshgrid(linspace(-m, m, N));

xdot = y;
ydot = -x + (1/3).* x.^3 - y;

%quiver(x, y, xdot, ydot)

%initial and unsafe sets
theta = linspace(0, 2*pi, 100);
circ = [cos(theta); sin(theta)];
X0 = [1.5; 0] + circ*0.5;
Xu = [-1; -1] + circ*0.4;
%p(x) = -(x+1)^2 + (y+1)^2

figure(3)
hold on
streamline(x, y, xdot, ydot, x, y)
plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3)
plot(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3)

hold off
axis square
xlabel('x')
ylabel('y')
title('Flow Field')