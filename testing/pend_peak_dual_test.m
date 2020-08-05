%pendulum with unit constants (length, mass, gravity)

%find maximum angular velocity with hybrid controller
%mset clear
rng(30, 'twister')

%states/variables
x = sdpvar(3, 1);
c = x(1);
s = x(2);
w = x(3);

%objective to maximize
objective = x(3)^2;        %maximum angular velocity
%objective = c+1;
%objective = 0.5*w^2 - c;    %maximum energy

%support
Xsupp.eq = (c^2 + s^2 - 1);

%% Formulate dynamics
%friction
b = 0.001;

energy = 0.5*w^2 - c;
energy_gap = energy - 1;
f_wait = [-s*w; c*w; -s - b*w];

%lqr realization
A = [0 1; 1 -b];
B = [0; -1];

Q = diag([2, 1]);
R = 1;
[K, S, clp] = lqr(A, B, Q, R);



%this is a stable choice when |th-pi|<=pi/4 and |w| <= 1.
u_lqr = @(x) K*[-sin(x(1)); x(2)];

u_lqr_poly = K*[-s; w];
%lqr decision
quad_lqr = [-s; w]'*S*[-s; w];

f_lqr = f_wait + [0; 0; 1]*u_lqr_poly;

%swing-up
k = 1; %swing-up gain
u_swing = @(x) -k * x(2) * ((0.5*x(2)^2 - cos(x(1))) - 1);

u_swing_poly = -k * w * energy_gap;
f_swing = f_wait + [0; 0; 1]*u_swing_poly;




%parameters of hybrid controller
epsilon = 0.1; %decision boundary for swing-up
delta = 2; %decision boundary for LQR-like control

%hybrid controller
u = @(x) u_hybrid(x, S, epsilon, delta, u_lqr, u_swing);


X_swing.ineq = (energy_gap^2 -  epsilon^2);
% X_swing_pos = energy_gap >= epsilon;
% X_swing_neg = energy_gap <= -epsilon;
X_wait.ineq  = [-energy_gap^2 + epsilon^2; quad_lqr - delta];
X_lqr.ineq   = [-energy_gap^2 + epsilon^2; -(quad_lqr - delta)];
% X_wait  = [energy_gap <= epsilon; energy_gap >= -epsilon; quad_lqr >= delta];
% X_lqr   = [energy_gap <= epsilon; energy_gap >= -epsilon; quad_lqr <= delta];


%wrap up the dynamics
% f = {f_swing, f_wait, f_lqr};
% X = {X_swing, X_wait, X_lqr};

%f = {f_lqr};%
%X = {[]};
%X = {(epsilon^2 - energy_gap^2)};

%this works. now to start changing the dynamics
%f = {f_lqr, f_lqr, f_lqr};
%X = {X_lqr, X_wait, X_swing};

f = {f_lqr, f_swing, f_swing};

%f = {f_lqr, f_wait, f_wait};
X = {X_lqr, X_wait, X_swing};


%initial state
%theta starts at pi, and goes to pi +- th_max
%th_max = pi;


%entirely in LQR portion
%th_max = pi/8;
%w_max = 0;

%th_max = pi/2;
%w_max = 0;
th_max = pi - 0.01;
w_max = 1.5;

HALF_ANGLE = 1; %keep sin(theta) >= 0

%X0 = [c <= cos(th_max); s^2 <= sin(th_max)^2; w^2 <= w_max^2];
%X0 = [c >= cos(th_max); s^2 <= sin(th_max)^2; w^2 <= w_max^2];
if HALF_ANGLE
    if w_max == 0
        X0.ineq = [-cos(th_max) - c; s];
        X0.eq = [c^2 + s^2 - 1, w];
    else
        X0.ineq = [-cos(th_max) - c; s; w_max^2 - w^2];
        X0.eq = [c^2 + s^2 - 1];        
    end
else
    if w_max == 0
        X0.ineq = [-cos(th_max) - c];
        X0.eq = [c^2 + s^2 - 1, w];
    else
        X0.ineq = [-cos(th_max) - c; w_max^2 - w^2];
        X0.eq = [c^2 + s^2 - 1];
    end    
end

%formulate problem
p_opt = peak_options;
p_opt.var.x = x;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.state_init = X0;
p_opt.state_supp = Xsupp;
options.R = 8;

%p_opt.scale = 0;

p_opt.obj = objective;

order = 7;
out = peak_estimate_dual(p_opt, order);
peak_val = out.peak_val;
ang_max = sqrt(-peak_val);

%% Simulate
Nsample = 100;
Tmax_sim = 20;

if HALF_ANGLE
    sampler = @() [pi - rand()* th_max; (2*rand()-1)*w_max];
else
    sampler = @() [pi + (2*rand() - 1)* th_max; (2*rand()-1)*w_max];
end
out_sim = pend_sampler(sampler, Nsample, Tmax_sim, out.dynamics.nonneg, 0, u);

if (out.optimal == 1)
    x0_angle = [atan2(out.x0(2), out.x0(1)); out.x0(3)];
    out_sim_peak = pend_sampler(x0_angle, 1, Tmax_sim, out.dynamics.nonneg, b);
    nplot = nonneg_plot(out_sim, out_sim_peak);
    pplot = pend_plot_3(out_sim, out_sim_peak);
%     splot = state_plot_3(out, out_sim, out_sim_peak);
% %     splot = state_plot_N(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out_sim);
%     splot = state_plot_2(out, out_sim);
%     splot3 = state_plot_3(out, out_sim);
    %splot = state_plot_N(out, out_sim);
    pplot = pend_plot_3(out_sim);
end

function u_out = u_hybrid(x, S, eps, delta, u_lqr, u_swing)
%hybrid controller.
%if abs(energy -1) >= eps, use u_swing
%if abs(energy-1) <= eps and norm(x)^2 <= delta, use u_lqr
%else, wait
energy_gap = 0.5*x(2)^2 - cos(x(1)) - 1;

if abs(energy_gap) >= eps
    u_out = u_swing(x);
else
    %if norm(x)^2 <= delta
    sx = [sin(x(1) - pi); x(2)];
    Sform = sx'*S*sx;
    if Sform <= delta
        u_out = u_lqr(x);
    else
        %u_out = 0;
        u_out = u_swing(x);
    end
end

end