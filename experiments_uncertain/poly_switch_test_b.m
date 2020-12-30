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
mpol('x', 2, 1);
mpol('w', 1, 1);
% mpol('d', 2, 1);
mpol('b', 2, 1);

%support
Xsupp = [];

%dynamics
%in principle, d is in a box
% dmax = 0.5;
% dmax = 0.4;
% dmax = 0.3;
% dmax = 0.2;
dmax = 2;

% wmax = 0.5;
draw = dmax.*(2.*b - 1);
w = 0.1;
f1 = [-(5 + draw(1))*x(1) - 8*x(2) + w*x(1)*x(2); - 2*x(2) - w*x(1)^3];
f2 = [-2*x(1) - 4*x(2) - w*x(2)^3; (10 + draw(2))*x(1) - 2*x(2) + w*x(1)*x(2)];

f = {f1, f2};
X = {[], []};


%initial set
C0 = -[1; 1];
R0 = 0.5;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%objective to maximize
objective = x(2);

%
p_opt = peak_options;
p_opt.var.t = t;
p_opt.var.x = x;
p_opt.var.b = b;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

Tmax_sim = 4;
if TIME_VARYING
    p_opt.Tmax = Tmax_sim;
end

% p_opt.box = 3;
p_opt.box = [-1.5, 2; -3, 3];
% p_opt.box = [-2.5, 2; -2.5, 2.5];
% p_opt.box = [-4, 0.5, 0; 2, 3.5, 4]';
p_opt.scale = 0;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

% order = 5;
order = 4;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;

% mu = 1;
Nsample = 500;

% Nsample = 150;
% Nsample = 50;
% Nsample = 20;
% sampler = @() circle_sample(1)'*R0 + C0;

s_opt = sampler_options;
s_opt.sample.x = @() sphere_sample(1, 2)'*R0 + C0;
s_opt.Tmax = Tmax_sim;
s_opt.Nb = 2;
s_opt.parallel = 1;


s_opt.mu = 0.2;
tic
out_sim = sampler(out.dynamics, Nsample, s_opt);
sample_time = toc;
disp(['Sampler Elapsed Time: ' , num2str(sample_time), ' seconds'])

    nplot = nonneg_plot(out, out_sim);

    splot = state_plot_2_lim(out, out_sim);
%     splot = state_plot_N(out, out_sim);
% end
