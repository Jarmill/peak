%test out arguments of peak_estimate routine
%clear

SOLVE = 1;
PLOT = 1;

if SOLVE

mset clear
% rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 2, 1);
mpol('d', 1, 1);

dmax = 0.2;

draw = dmax*(2*d - 1);

%support
Xsupp = [];

%dynamics


% A = [-0.3095    0.8900;
%    -0.7408   -0.3253];
A = [-0.3 0.8;
    -0.75 -0.3];
f1 = [A*x] + [0; draw];
X1 = [];

f2 = [0.8 0.5; -0.5 0.8]*x + [0; draw];
X2 = [x(1) >= 0];

POLY = 1;
if POLY
%     f1 = f1 + 0.1*[-x(1)*x(2); x(1)*x(2)];
    f2 = f2 + 0.01*[-x(1)^2; -x(1)*x(2)];
end


% f = f2;
% X = X2;

f = {f1, f2};
X = {X1, X2};


%initial set
%C0 = [1.2; 0];
C0 = [-1.5; 0];
R0 = 0.4;

%C0 = [0.8; 0];
%R0 = 0.1;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%objective to maximize
objective = -x(2);
%objective = -x(2) - x(1);

% objective = [-x(2); -x(1); -x(1) - x(2)];
% objective = [-x(2); -x(1)];
%
p_opt = peak_options;
p_opt.var.x = x;
p_opt.var.d = d;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;
p_opt.disturb = (d*(1-d) >= 0);

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.dynamics.discrete = 1;

Tmax_sim = 30;
p_opt.Tmax = Tmax_sim;

p_opt.box = 3;
% p_opt.box = [-2, 2; -2, 1.5];
p_opt.scale = 0;
%p_opt.R = 6;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 4;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;

end
if PLOT
%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;

mu = 1;

Nsample = 100;
Nsample = 250;
s_opt = sampler_options;
s_opt.sample.x = @() sphere_sample(1, 2)'*R0 + C0;
s_opt.sample.d = @() dmax * (2*rand() - 1);
s_opt.Nd = 1;
s_opt.Tmax = Tmax_sim;
s_opt.parallel = 0;
out_sim = sampler(out.dynamics, Nsample, s_opt);

% if (out.optimal == 1)
%     s_opt.sample.x = out.x0;
%     out_sim_peak = sampler(out.dynamics, 1, s_opt);
%     nplot = nonneg_plot(out, out_sim, out_sim_peak);
%     cplot = cost_plot(out, out_sim, out_sim_peak);
%     splot = state_plot_2(out, out_sim, out_sim_peak);
%     %     splot = state_plot_N(out, out_sim, out_sim_peak);
% else
    nplot = nonneg_plot(out, out_sim);
    cplot = cost_plot(out, out_sim);
%     splot = state_plot_2(out, out_sim);
    splot = state_plot_2_lim(out, out_sim);
%     splot = state_plot_N(out, out_sim);
% end
end