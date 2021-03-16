%test out arguments of peak_estimate routine
%clear
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 1, 1);

%support
% Xsupp = [x >= 0];
% Xsupp = (x^2 <= 1);
Xsupp = [];
%dynamics
% f = [-x(1)*(1-x(3)); -x(2)*x(3); -x(3)];
f = x^2 - x^3;
X = [];

%initial set
% C0 = [1.5; 0.5; 0.5];
% R0 = 0.4;


% X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);
% X0 = -x*(x+1) >= 0;
X0 = -(x+0.1)*(x+1) >= 0;
% X0 = (x + 1)
%objective to maximize
%objective = -x(2);
objective = 4*x*(1-x);


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

p_opt.box = 1;
p_opt.scale = 0;
%p_opt.R = 6;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 5;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

% Nsample = 100;
% Nsample = 20;
% sampler = @() ball_sample(1, 3)'*R0 + C0;
% 
% out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, 0, @ode45);
% 
% if (out.optimal == 1)
%     out_sim_peak = switch_sampler(out.dynamics, out.x0, 1, Tmax_sim);
%     nplot = nonneg_plot(out, out_sim, out_sim_peak);
% 
%     splot = state_plot_3(out, out_sim, out_sim_peak);
%     splot = state_plot_N(out, out_sim, out_sim_peak);
% else
%     nplot = nonneg_plot(out, out_sim);
%     splot = state_plot_3(out, out_sim);
%     splot = state_plot_N(out, out_sim);
% end
