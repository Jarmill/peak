%test out arguments of peak_estimate routine
%phytoplankton growth model
%example 5.5 of https://arxiv.org/pdf/1703.05085.pdf

SOLVE = 1;
PLOT = 1;

if SOLVE

mset clear
% rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 3, 1);


%dynamics
f = x + 0.01*[  1-x(1)-0.25*x(1)*x(2); 
                x(2)*(2*x(3)-1); 
                0.25*x(1)-2*x(3)^2];
X = [];


%initial set
% X0 = [(x(1) - 0.25)^2 <= 0.05^2; (x(2) - 0.25)^2 <= 0.05^2; x(3)^2 <= 0.05^2];

X0 = [x(1) >= -0.3; x(1) <= -0.2; x(2) >= -0.3; x(2) <= -0.2; x(3) >= -0.05; x(3) <= 0.05];

%objective to maximize
objective = -x(3);

p_opt = peak_options;
p_opt.var.x = x;
% p_opt.var.d = d;

p_opt.state_init = X0;


p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.dynamics.discrete = 1;

Tmax_sim = 50;
p_opt.Tmax = Tmax_sim;

% p_opt.box = 3;
% p_opt.box = [-2, 2; -2, 1.5];

p_opt.box = [-0.5, 1.5; -0.5, 0.5; -0.5, 0.5];
p_opt.scale = 0;
%p_opt.R = 6;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 3;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;

end
if PLOT
%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
mu = 1;

Nsample = 100;
% Nsample = 250;
s_opt = sampler_options;
s_opt.sample.x = @() rand(3, 1).*[0.1; 0.1; 0.05] + [-0.3; 0.3; -0.05];
% s_opt.Nd = 1;
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
    splot = state_plot_3(out, out_sim);
%     splot = state_plot_2_lim(out, out_sim);
%     splot = state_plot_N(out, out_sim);
% end
end