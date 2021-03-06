%test out arguments of peak_estimate routine
%clear
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

TIME_VARYING = 1;

%define variables
%gloptipoly is sensitive to the order in which variables are declared when
%defining a measure

mpol('t', 1, 1);
mpol('x', 3, 1);
mpol('b', 2, 1);
mpol('w', 1, 1);

%support
Xsupp = [];

%dynamics
%in principle, d is in a box
% dmax = 0.5;
% dmax = 0.4;
% dmax = 0.3;
% dmax = 0.2;
% dmax = 0.1;
% dmax=0;
% draw = dmax.*(2.*d - 1);

wmax = 0.1;

dmax = [0.5; 0.2];

% d = dmax *(2*b - 1);
d = [dmax.*(2*b-1)];

a0 = 1 + d(1);
b0 = 0.5 + d(2);
% G0 = 2 + wmax*w;
G0 = 2;
f0 = [a0*x(1) + b0*x(2) + x(3) - 2*x(2)^2;
    a0*x(2) - b0*x(1) + 2*x(1)*x(2);
    -G0*x(3) - 2*x(1)*x(3)];

f00 = subs(f0, b, [0; 0]);
f10 = subs(f0, b, [1; 0]);
f01 = subs(f0, b, [0; 1]);
f11 = subs(f0, b, [1; 1]);

% X = [];
f = {f00, f10, f01, f11};
X = {[],[],[],[]};


%initial set
C0 = [1; 1; 1];
R0 = 0.4;

%C0 = [0.8; 0];
%R0 = 0.1;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 + (x(3)-C0(3))^2 <= R0^2);
%objective to maximize
objective = x(2);
%objective = -x(2) - x(1);
%
p_opt = peak_options;
p_opt.var.t = t;
p_opt.var.x = x;
% p_opt.var.w = w;
% p_opt.var.b = b;

% p_opt.var.b = b;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;
% p_opt.disturb = (d.^2 <= dmax^2);
% p_opt.param = (w^2 <= 1);

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

Tmax_sim = 5;
if TIME_VARYING
    p_opt.Tmax = Tmax_sim;
end

% p_opt.box = 4;
% p_opt.box = [-1, 3; -1.5, 2];
p_opt.box = [-4, 0.5, 0; 3, 3.6, 4]';
p_opt.scale = 0;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 3;
% order = 2;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;

% mu = 1;

Nsample = 250;
% Nsample = 50;
% Nsample = 20;
% sampler = @() circle_sample(1)'*R0 + C0;

s_opt = sampler_options;
% s_opt.sample.x = @() circle_sample(1)'*R0 + C0;
s_opt.sample.x = @() sphere_sample(1, 3)'*R0 + C0;
% s_opt.sample.w = @() wmax * (2*rand(1, 1) - 1);
s_opt.Tmax = Tmax_sim;
% s_opt.Nb = 2;
% s_opt.Nw = 1;
s_opt.parallel = 1;


s_opt.mu = 0.2;
tic
out_sim = sampler(out.dynamics, Nsample, s_opt);
sample_time = toc;
disp(['Sampler Elapsed Time: ' , num2str(sample_time), ' seconds'])
% out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, 0, @ode45);

% if (out.optimal == 1)
%     out_sim_peak = switch_sampler(out.dynamics, out.x0, 1, Tmax_sim);
%     nplot = nonneg_plot(out_sim, out_sim_peak);
% 
%     splot = state_plot_2(out, out_sim, out_sim_peak);
%     %     splot = state_plot_N(out, out_sim, out_sim_peak);
% else
    nplot = nonneg_plot(out, out_sim);
    splot = state_plot_3(out, out_sim);
%     splot = state_plot_N(out, out_sim);
% end
