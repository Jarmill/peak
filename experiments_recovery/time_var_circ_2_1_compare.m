%Example 2.1 from the Fantuzzi/Goluskin paper
%on a circle
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 2, 1);
mpol('t', 1, 1);
%support
Xsupp = [];

%dynamics
f = [x(2)*t - 0.1*x(1) - x(1)*x(2);
    -x(1)*t - x(2) + x(1)^2];
X = [];

%initial set
C0 = [-0.75; 0];
R0 = 1;

EQ = 1;

if EQ
    X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 == R0^2);
    sampler = @() sphere_sample(1, 2)'*R0 + C0;
else
    X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);
    sampler = @() ball_sample(1, 2)'*R0 + C0;
end
%objective to maximize
objective = x(1);
%objective = [x(1);  1/sqrt(5) * x(1) + 2/sqrt(5) *x(2)];

% objective = [x(1); x(2)];
% 
p_opt = peak_options;
p_opt.var.x = x;
p_opt.var.t = t;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

Tmax_sim = 5;
p_opt.Tmax = Tmax_sim;

p_opt.box = [-3, 2; -2, 2];
p_opt.scale = 0;
%p_opt.R = 6;


p_opt.rank_tol = 1e-3;
p_opt.obj = objective;

% order = 3;
% out = peak_estimate(p_opt, order);
% peak_val = out.peak_val;

degree_bound = (5);
peak_fin_horizon = zeros(length(degree_bound), 1);

out_fin_horizon = cell(length(degree_bound), 1);


%run experiments
for i = 1:length(degree_bound)
    order = degree_bound(i);
    out_fin_horizon{i} = peak_estimate(p_opt, order);
    peak_fin_horizon(i) = out_fin_horizon{i}.peak_val;
end

