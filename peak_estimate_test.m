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
f1 = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
f2 = [-x(1); -x(2)];
f3 = [-0.1*x(1)-x(2); x(1) - 0.1*x(2)];
%X1 = [x(2) <= 2];
X1 = [];
X2 = [];
X3 = [];
%X3 = [x(1) >= 0];

%f = {f1, f2, f3};
%X = {X1, X2, X3};

f = f1;
X = X1;

%f = {f1, f3};
%X = {X1, X3};
%


%f = {f3, f2};
%X = {X3, X2};


%f4 = [0 2; -1 -1]*x;
%f5 = [1 2; -3 -2]*x;
%f = {f4, f5};
%X = {[], []};

%f = f1;
%X = X1;
% 
% f = f3;
% X = X3;

if iscell(f)
    nsys = length(f);
else
    nsys = 1;
end


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

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

Tmax_sim = 5;
%p_opt.Tmax = Tmax_sim;
p_opt.state_init = X0;
p_opt.state_supp = Xsupp;
%p_opt.box = [-4, 6];
p_opt.box = 3;
p_opt.scale = 0;
%p_opt.R = 6;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 4;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;

mu = 1;

%Nsample = 100;
Nsample = 10;
sampler = @() circle_sample(1)'*R0 + C0;



out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, @ode45);

if (out.optimal == 1)
    out_sim_peak = switch_sampler(out.dynamics, out.x0, 1, Tmax_sim);
    nplot = nonneg_plot(out_sim, out_sim_peak);

    splot = state_plot_2(out, out_sim, out_sim_peak);
    %     splot = state_plot_N(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out_sim);
    splot = state_plot_2(out, out_sim);
%     splot = state_plot_N(out, out_sim);
end
