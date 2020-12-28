%test out arguments of peak_estimate routine
%clear
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol


mpol('x', 2, 1);
mpol('d', 1, 1);

%support
Xsupp = [];

%dynamics
%in principle, d is in a box [0,1]
dmax = 0.4;
draw = dmax*(2*d - 1);

f = [x(2); -x(1) - x(2) + (1/3).* x(1).^3 + draw];
X = [];



%initial set
%C0 = [1.2; 0];
C0 = [1.5; 0];
R0 = 0.4;

%C0 = [0.8; 0];
%R0 = 0.1;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%objective to maximize
objective = -x(2);
%objective = -x(2) - x(1);
%
p_opt = peak_options;
p_opt.var.x = x;
p_opt.var.d = d;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;
% p_opt.disturb = (d^2 <= dmax^2);
p_opt.disturb = (d*(1-d) >= 0);

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;


Tmax_sim = 5;
% p_opt.Tmax = Tmax_sim;



% p_opt.box = 3;
p_opt.box = [-1, 3; -1.5, 2];
p_opt.scale = 0;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 5;
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
s_opt.Nd = 1;
s_opt.Tmax = Tmax_sim;

s_opt.parallel = 1;
s_opt.mu = 0.4;

out_sim = sampler(out.dynamics, Nsample, s_opt);

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
