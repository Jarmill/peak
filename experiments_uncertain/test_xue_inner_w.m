%test out arguments of peak_estimate routine
%clear
mset clear
rng(300, 'twister')

%dynamics and support set
%modified inner RoA Xue flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol

mpol('t', 1, 1);
mpol('x', 2, 1);
mpol('d', 1, 1);
mpol('w', 1, 1);

%order 4
%obj = 0.768

%support
Xsupp = [];

%dynamics
%in principle, d is in a box [0,1]
dmax = 0.2;
wmax = 0.5;
% draw = dmax*(2*d - 1);

% f = [x(2); -x(1) - x(2) + (1/3).* x(1).^3 + draw];

f = [-0.5*x(1) - (0.5 + d)*x(2) + 0.5; -0.5*x(2) + 1 + w];
X = [];



%initial set
%C0 = [1.2; 0];
C0 = [-1; -1];
R0 = 0.5;

%C0 = [0.8; 0];
%R0 = 0.1;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%objective to maximize
objective = x(1);

p_opt = peak_options;
p_opt.var.t = t;
p_opt.var.x = x;
p_opt.var.d = d;
p_opt.var.w = w;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;
p_opt.disturb = (d^2 <= dmax^2);
p_opt.param = (w^2 <= wmax^2);

% p_opt.disturb = (d*(1-d) >= 0);


p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;


Tmax_sim = 10;
p_opt.Tmax = Tmax_sim;



% p_opt.box = 3;
p_opt.box = [-2, 2; -2, 4];
p_opt.scale = 0;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 4;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;

% mu = 1;

Nsample = 100;
% Nsample = 50;
% sampler = @() circle_sample(1)'*R0 + C0;

s_opt = sampler_options;
s_opt.sample.x = @() sphere_sample(1, 2)'*R0 + C0;
s_opt.sample.d = @() dmax * (2*rand() - 1);
s_opt.sample.w = @() wmax * (2*rand() - 1);
s_opt.Nd = 1;
s_opt.Tmax = Tmax_sim;

s_opt.parallel = 0;
s_opt.mu = 0.4;

out_sim = sampler(out.dynamics, Nsample, s_opt);
out.optimal = 0;

% out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, 0, @ode45);

if (out.optimal == 1)
    s_opt.sample.x = out.x0;
    out_sim_peak = sampler(out.dynamics, 1, s_opt);
%     out_sim_peak = switch_sampler(out.dynamics, out.x0, 1, Tmax_sim);
    nplot = nonneg_plot(out, out_sim, out_sim_peak);
% 
    splot = state_plot_2(out, out_sim, out_sim_peak);
    %     splot = state_plot_N(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
    splot = state_plot_2(out, out_sim);
%     splot = state_plot_N(out, out_sim);
end
