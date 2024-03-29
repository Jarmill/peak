%test out arguments of peak_estimate routine
%clear
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
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
X = [];

%initial set
%C0 = [1.2; 0];
C0 = [1.5; 0];
R0 = 0.4;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%objective to maximize
objective = -x(2);
%objective = -x(2) - x(1);

% objective = [-x(2); -x(1); -x(1) - x(2)];
% objective = [-x(2); -x(1)];
%
p_opt = peak_options;
p_opt.var.x = x;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;
p_opt.dynamics.discrete = 0;

Tmax_sim = 10;
p_opt.Tmax = Tmax_sim;

p_opt.box = 4;
p_opt.scale = 0;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 1;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;

mu = 1;

Nsample = 100;
% Nsample = 20;
sampler = @() circle_sample(1)'*R0 + C0;

dynamics = struct('f', {@(x_in) eval(f, x, x_in)}, 'event', {@(t, x) all_event(t, x)});
out_sim = switch_sampler(dynamics, sampler, Nsample, Tmax_sim, 1, 0, @ode45);

if (out.optimal == 1)
    out_sim_peak = switch_sampler(dynamics, out.x0, 1, Tmax);
    nplot = nonneg_plot(out, out_sim, out_sim_peak);

    splot = state_plot_2(out, out_sim, out_sim_peak);
    %     splot = state_plot_N(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
    splot = state_plot_2(out, out_sim);
%     splot = state_plot_N(out, out_sim);
end


