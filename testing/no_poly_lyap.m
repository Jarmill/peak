%Example 2.1 from the Fantuzzi/Goluskin paper
%on a circle
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('t', 1, 1);
mpol('x', 2, 1);
%support
Xsupp = [];

%dynamics
f0 = [x(1)*(x(2)-1); -x(2)];

B = 10;
Q = diag([1,1]./B); % Scale to unit box
f = Q*subs(f0, x, inv(Q)*x);

X = [];

%initial set
C0 = [2; 1.5]/B;
R0 = 0.5/B;

EQ = 1;

if EQ
    X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 == R0^2);
    sampler_x = @() sphere_sample(1, 2)'*R0 + C0;
else
    X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);
    sampler_x = @() ball_sample(1, 2)'*R0 + C0;
end
%objective to maximize
objective = (B*x(1))^2;
%objective = [x(1);  1/sqrt(5) * x(1) + 2/sqrt(5) *x(2)];

% objective = [x(1); x(2)];
% 
p_opt = peak_options;
p_opt.var.x = x;
% p_opt.var.t = t;

%p_opt.state_supp = Xsupp;
p_opt.state_init = X0;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;
Tmax_sim = 5;
% p_opt.Tmax = Tmax_sim;
p_opt.Tmax = Inf;

p_opt.box = [-2, 2; -2, 2];


p_opt.rank_tol = 1e-3;
p_opt.obj = objective;

order = 3;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;

mu = 1;

% Nsample = 100;

Nsample = 50;

% Nsample = 20;


s_opt = sampler_options;
s_opt.sample.x = sampler_x;
s_opt.Tmax = Tmax_sim;

out_sim = sampler(out.dynamics, Nsample, s_opt);

% out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, 0, @ode45);

if (out.optimal == 1)
%     out_sim_peak = switch_sampler(out.dynamics, out.x0, 1, Tmax_sim);
    s_opt.sample.x = out.x0;
    out_sim_peak = sampler(out.dynamics, 1, s_opt);
    nplot = nonneg_plot(out, out_sim, out_sim_peak);
    cplot = cost_plot(out, out_sim, out_sim_peak);
    splot = state_plot_2(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
    cplot = cost_plot(out, out_sim);
    splot = state_plot_2(out, out_sim);
end
