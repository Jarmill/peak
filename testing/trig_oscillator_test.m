%pendulum with unit constants (length, mass, gravity)

mpol('x', 3, 1);
%states
c = x(1);   %cos(theta)
s = x(2);   %sin(theta)
w = x(3);   %angular velocity


%starting angular velocity

%support space
Xsupp = [c^2 + s^2 == 1];
%Xsupp = [c^2 == 1 - s^2];

%dynamics
f = [-s*w; c*w; 0];
X = [];


th_max = pi/3;
w_max = 3;
%X0 = [c <= cos(th_max); s^2 <= sin(th_max)^2; w^2 <= w_max^2];
X0 = [c >= cos(th_max); s^2 <= sin(th_max)^2; w^2 <= w_max^2];
% X0 = [c >= cos(th_max); s^2 <= sin(th_max)^2; w == w_start];

%maximum height of oscillator
objective = s;

%formulate problem
p_opt = peak_options;
p_opt.var.x = x;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.state_init = X0;
p_opt.state_supp = Xsupp;
p_opt.box = [-1, 1; -1, 1; -4, 4];
p_opt.scale = 0;

p_opt.obj = objective;

%run problem 
order = 2;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;



Tmax_sim = 10;
Nsample = 20;
%sampler = @() circle_sample(1)'*R0 + C0;

%sampler = @() pend_sampler(th_max, 0) + [0; 0; w_start];
sampler = @() pend_sampler(th_max, w_max);

mu = 1;
out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, @ode45);
% 
if (out.optimal == 1)
    out_sim_peak = switch_sampler(out.dynamics, out.x0, 1, Tmax_sim);
%     nplot = nonneg_plot(out_sim, out_sim_peak);
% 
%     splot = state_plot_2(out, out_sim, out_sim_peak);
        splot = state_plot_N(out, out_sim, out_sim_peak);
else
%     nplot = nonneg_plot(out_sim);
%     splot = state_plot_2(out, out_sim);
    splot = state_plot_N(out, out_sim);
end

