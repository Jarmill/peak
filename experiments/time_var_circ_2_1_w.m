%Example 2.1 from the Fantuzzi/Goluskin paper
%on a circle
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

PARAM = 1;

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 2, 1);
mpol('w', 1, 1);
mpol('t', 1, 1);

%support
Xsupp = [];

%dynamics
f = [x(2)*t - (0.1 + w)*x(1) - x(1)*x(2);
    -x(1)*t - x(2) + x(1)^2];
X = [];


% w_max = 0.05;
% Wsupp = [w^2 <= w_max^2];
w_max = 0.05;
Wsupp = [w^2 <= w_max^2];
%initial set
C0 = [-0.75; 0];
R0 = 1;

EQ = 1;

if EQ
    X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 == R0^2);    
    sampler = @() [sphere_sample(1, 2)'*R0 + C0; (2*rand(1) - 1)*w_max];
else
    X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);
    sampler = @() [ball_sample(1, 2)'*R0 + C0;(2*rand(1) - 1)*w_max];
end


%objective to maximize
objective = x(1);
%objective = [x(1);  1/sqrt(5) * x(1) + 2/sqrt(5) *x(2)];


%there's a bug in beta'*cost. cost should be cost_all. Fix this.
% objective = [x(1); x(2)];

p_opt = peak_options;
p_opt.var.x = x;
p_opt.var.t = t;
p_opt.var.w = w;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;
p_opt.param = Wsupp;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;
p_opt.dynamics.discrete = 0;
Tmax_sim = 5;
p_opt.Tmax = Tmax_sim;
% 
% p_opt.box = [-3, 2; -2, 2];

p_opt.box = 3;
p_opt.scale = 0;
%p_opt.R = 6;


p_opt.rank_tol = 1e-3;
p_opt.obj = objective;

order = 3;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;

mu = 1;

Nsample = 100;

% Nsample = 50;

% Nsample = 20;

out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, 1, @ode45);

if (out.optimal == 1)
    out_sim_peak = switch_sampler(out.dynamics, [out.x0; out.w], 1, Tmax_sim, mu, 1);
    nplot = nonneg_plot(out, out_sim, out_sim_peak);
    cplot = cost_plot(out, out_sim, out_sim_peak);
    splot = state_plot_2(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
    cplot = cost_plot(out, out_sim);
    splot = state_plot_2(out, out_sim);
end
