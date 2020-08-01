%test out arguments of peak_estimate routine
%clear
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
%mpol('x', 2, 1);
x = sdpvar(2, 1);

%support
Xsupp = [];

%dynamics
f1 = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)]; %prajna and rantzer
f2 = [-x(1); -x(2)]; %go to origin
f3 = [-0.1*x(1)-x(2); x(1) - 0.1*x(2)]; %decaying spiral
X1 = [];
X2 = [];
X3 = [];
%X3 = [x(1) >= 0];

% f = {f1, f2, f3};
% X = {X1, X2, X3};
f = {f1, f2};
X = {X1, X2};

% f = f1;
% X = X1;

f4 = [0 2; -1 -1]*x;
f5 = [1 2; -3 -2]*x;
f = {f4, f5};
X = {[], []};

%initial set

%R0 = 0.4;
%C0 = [1.5; 0];%

R0 = 0.5;
C0 = [-1; -1];

X0 = struct;
%same trick, + instead of - in X0
%that's where the other bug was.
X0.ineq = R0^2 - (x(1)-C0(1))^2 - (x(2)-C0(2))^2;

%objective to maximize
objective = -x(2);
%objective = -x(2) - x(1);
%
p_opt = peak_options;
p_opt.var.x = x;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

Tmax_sim = 10;
%p_opt.Tmax = Tmax_sim;
p_opt.state_init = X0;
p_opt.state_supp = Xsupp;
p_opt.scale = 0;
p_opt.R = 5;

p_opt.obj = objective;

order = 5;
out = peak_estimate_dual(p_opt, order);
peak_val = out.peak_val;

% if ~out.feasible
%     p_opt.prev_cost = peak_val;
%     out2 = peak_estimate_dual(p_opt, order);
%     peak_val_2 = out2.peak_val;
% end


%% now do plots
%there is something wrong with the definition of v. why?


rng(50, 'twister')
x0 = C0;

mu = 1;

%Nsample = 100;
Nsample = 50;
sampler = @() circle_sample(1)'*R0 + C0;
% 
% 
% 
out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu);
% 
% if (out.optimal == 1)
%     out_sim_peak = switch_sampler(out.dynamics, out.x0, 1, Tmax_sim);
%     nplot = nonneg_plot(out_sim, out_sim_peak);
% 
% %     splot = state_plot_2(out, out_sim, out_sim_peak);
%         splot = state_plot_N(out, out_sim, out_sim_peak);
% else
nplot = nonneg_plot(out_sim);
splot = state_plot_2(out, out_sim);
% splot = state_plot_N(out, out_sim);
% end
