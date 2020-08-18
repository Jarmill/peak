%Example 2.1 from the Fantuzzi/Goluskin paper
%on a point
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
C0 = [0; 1];
sampler = C0;
X0 = (x==C0);
%objective to maximize
objective = x(1);
%objective = [x(1);  1/sqrt(5) * x(1) + 2/sqrt(5) *x(2)];

%objective = [x(1); x(2)];

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

order = 3;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')

mu = 1;

% Nsample = 100;

% Nsample = 50;
Nsample = 1;

% Nsample = 20;

out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, 0, @ode45);

if (out.optimal == 1)
    out_sim_peak = switch_sampler(out.dynamics, out.x0, 1, Tmax_sim);
    nplot = nonneg_plot(out, out_sim, out_sim_peak);
    cplot = cost_plot(out, out_sim, out_sim_peak);
    splot = state_plot_2(out,  out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
    cplot = cost_plot(out, out_sim);
    splot = state_plot_2(out, out_sim);
end
