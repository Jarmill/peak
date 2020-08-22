%pendulum with unit constants (length, mass, gravity)

mpol('x', 3, 1);
%states
c = x(1);   %cos(theta)
s = x(2);   %sin(theta)
w = x(3);   %angular velocity

%a little friction
%b = 0.0;
% b = 0.1;
b = 0.05;

%open loop control for now
u = 0; 

%support space
% Xsupp = (c^2 + s^2 <= 1); %easier problem (less space to be >= 0), but
%v(x) is inaccurate

Xsupp = [c^2 + s^2 == 1];


% Xsupp = []; %harder problem, but accurate v(x)

%dynamics
f = [-s*w; c*w; -s - b*w];
X = [];
% X = (c^2 + s^2 == 1);


% th_max = pi/4;
th_max = pi/2;
% w_max = 0.0;
w_max = 1;

HALF_ANGLE = 0; %keep sin(theta) >= 0

%X0 = [c <= cos(th_max); s^2 <= sin(th_max)^2; w^2 <= w_max^2];
%X0 = [c >= cos(th_max); s^2 <= sin(th_max)^2; w^2 <= w_max^2];
if HALF_ANGLE
%     X0 = [c >= cos(th_max); s <= sin(th_max); s >= 0; w^2 <= w_max^2];
    X0 = [c >= cos(th_max); s >= 0; w^2 <= w_max^2];
else
    X0 = [c >= cos(th_max); w^2 <= w_max^2];
end
% X0 = [c >= cos(th_max); s^2 <= sin(th_max)^2; w^2 <= w_max^2; c^2 + s^2 == 1];
% X0 = [X0];

%maximum height of pendulum
%objective = 1-c;
objective = -c;
Xsupp = [Xsupp; (w^2 <= 9)];

%formulate problem
p_opt = peak_options;
p_opt.var.x = x;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.state_init = X0;
p_opt.state_supp = Xsupp;
p_opt.box = [];
p_opt.scale = 0;

p_opt.obj = objective;

%run problem 
order = 4;
out = peak_estimate(p_opt, order);
pend_height = 1-out.peak_val; %min vs. max?



Tmax_sim = 10;
Nsample = 20;
if HALF_ANGLE
    sampler = @() [rand(1, 1)*th_max; (2*rand(1, 1) - 1)*w_max];
else
    sampler = @() (2*rand(2, 1) - 1).*[th_max; w_max];
end

out_sim = pend_sampler(sampler, Nsample, Tmax_sim, out.dynamics.nonneg, b);


% % 
if (out.optimal == 1)
    x0_angle = [atan2(out.x0(2), out.x0(1)); out.x0(3)];
    out_sim_peak = pend_sampler(x0_angle, 1, Tmax_sim, out.dynamics.nonneg, b);
    nplot = nonneg_plot(out, out_sim, out_sim_peak);
% 
    splot = state_plot_3(out, out_sim, out_sim_peak);
% %     splot = state_plot_N(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
%     splot = state_plot_2(out, out_sim);
    splot3 = state_plot_3(out, out_sim);
%     splot = state_plot_N(out, out_sim);
end
% 
