%test out arguments of peak_estimate routine
%clear
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 3, 1);

sigma = 10;
rho = 28;
beta = 8/3;

%support
Xsupp = [];

%dynamics
%f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

f = [sigma*(x(2) - x(1));
    x(1)*(rho - x(3) - x(2));
    x(1)*x(2) - beta*x(3)];
X = [];

%initial set
%C0 = [1.2; 0];
C0 = [1; 1; 1];
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

Tmax_sim = 30;
p_opt.Tmax = Tmax_sim;

p_opt.box = 15;
p_opt.scale = 0;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 2;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;

mu = 1;

% Nsample = 100;
Nsample = 20;
sampler = @() sphere_sample(1, 3)'*R0 + C0;

out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, 0, @ode45);

if (out.optimal == 1)
    out_sim_peak = switch_sampler(out.dynamics, out.x0, 1, Tmax_sim);
    nplot = nonneg_plot(out, out_sim, out_sim_peak);

    splot = state_plot_3(out, out_sim, out_sim_peak);
    %     splot = state_plot_N(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
    splot = state_plot_3(out, out_sim);
%     splot = state_plot_N(out, out_sim);
end