%Example 2.1 from the Fantuzzi/Goluskin paper
%on a circle
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 2, 1);

%support
Xsupp = [];

%dynamics
% f = [x(2)*t - 0.1*x(1) - x(1)*x(2);
%     -x(1)*t - x(2) + x(1)^2];

r2 = x'*x;

A = [0.2, 1; 0, -0.4];
J = [0, -1; 1,  0];
f = A*x + J*x*r2;

X = [];

% X = (x(2) >= 0)
%initial set
C0 = [0; 0];
R0 = 0.5;

EQ = 1;
% SYM = 0;

if EQ
    X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 == R0^2);
    sampler = @() sphere_sample(1, 2)'*R0 + C0;
else
    X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);
    sampler = @() ball_sample(1, 2)'*R0 + C0;
end

% if ~SYM
%     sampler0 = sampler;
%     X = [X; V_sym'*x >= 0];
%     sampler = sample_rejector(V_sym, sampler0);
% end


%objective to maximize
% objective = x(1);
%objective = [x(1);  1/sqrt(5) * x(1) + 2/sqrt(5) *x(2)];

objective = r2;

% objective = [x(1); x(2)];

p_opt = peak_options;
p_opt.var.x = x;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

Tmax_sim = 20;
% p_opt.Tmax = Tmax_sim;

p_opt.box = 2;
p_opt.scale = 0;
%p_opt.R = 6;


p_opt.rank_tol = 3e-4;
p_opt.obj = objective;

% order = 7;
% out = peak_estimate(p_opt, order);
% peak_val = out.peak_val;    

degree_bound = (2:9);
peak_inf_horizon = zeros(length(degree_bound), 1);
peak_fin_horizon = zeros(length(degree_bound), 1);

out_inf_horizon = cell(length(degree_bound), 1);
out_fin_horizon = cell(length(degree_bound), 1);


%run experiments

p_opt_fin = p_opt;
p_opt_fin.Tmax = 3;

for i = 1:length(degree_bound)
    order = degree_bound(i);
    out_inf_horizon{i} = peak_estimate(p_opt, order);
    peak_inf_horizon(i) = out_inf_horizon{i}.peak_val;

    out_fin_horizon{i} = peak_estimate(p_opt, order);
    peak_fin_horizon(i) = out_fin_horizon{i}.peak_val;
end

save('sym_attractor_experiment', 'degree_bound', 'peak_inf_horizon', 'peak_fin_horizon')