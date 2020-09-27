%Example 2.1 from the Fantuzzi/Goluskin paper
%on a circle
mset clear
rng(300, 'twister')

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
A = [0.2, 1; 0, -0.4];
J = [0, -1; 1,  0];
f = @(t, x) A*x + J*x*(x'*x);

C0 = [0; 0];
R0 = 0.5;

EQ = 1;

Nsample = 24;
theta_sample = ((1:Nsample)-1)*2*pi/(Nsample);
% theta_sample = linspace(0, 2*pi, 16);

theta = linspace(0, 2*pi, 200);

figure(1)
clf

Tmax_sim = 20;
hold on
for i = 1:Nsample
    px = R0*[cos(theta_sample(i)), sin(theta_sample(i))];
    [T, X] = ode45(f, [0, Tmax_sim], px');
    
    plot(X(:, 1), X(:, 2), 'c')
end
plot(R0*cos(theta), R0*sin(theta), 'k', 'Linewidth', 2)
scatter(R0*cos(theta_sample), R0*sin(theta_sample), 50,  'ok')
xlabel('x_1')
ylabel('x_2')
xlim([-1.5,1.5])
ylim([-1.5,1.5])
axis square
hold off