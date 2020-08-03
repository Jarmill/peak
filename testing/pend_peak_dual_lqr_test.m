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
%no friction

%lqr realization
A = [0 1; 1 0];
B = [0; -1];

Q = diag([2, 1]);
R = 1;
[K, S, clp] = lqr(A, B, Q, R);



%this is a stable choice when |th-pi|<=pi/4 and |w| <= 1.
u_lqr = @(x) K*[-sin(x(1)); x(2)];

u_lqr_poly = K*[-s; w];
%lqr decision
quad_lqr = [-s; w]'*S*[-s; w];


f_wait = [-s*w; c*w; -s];
f_lqr = f_wait + [0; 0; 1]*u_lqr_poly;




%wrap up the dynamics
%f = {f_swing, f_wait, f_lqr};
%X = {X_swing, X_wait, X_lqr};
f = f_lqr;
X = [];

%initial state
%theta starts at pi, and goes to pi +- th_max
%th_max = pi;
%th_max = pi/4;
%w_max = 1;

th_max = pi/8;
w_max = 0;

HALF_ANGLE = 0; %keep sin(theta) >= 0

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
p_opt.R = 5;
%p_opt.Tmax = 10;

p_opt.obj = objective;

order = 4;
out = peak_estimate_dual(p_opt, order);
peak_val = out.peak_val;
ang_max = sqrt(-peak_val);

%% Simulate
Nsample = 20;
Tmax_sim = 10;

sampler = @() [pi + (2*rand() - 1)* th_max; (2*rand()-1)*w_max];

u = u_lqr;
out_sim = pend_sampler(sampler, Nsample, Tmax_sim, out.dynamics.nonneg, 0, u);

if (out.optimal == 1)
    x0_angle = [atan2(out.x0(2), out.x0(1)); out.x0(3)];
    out_sim_peak = pend_sampler(x0_angle, 1, Tmax_sim, out.dynamics.nonneg, b);
    splot = state_plot_N(out, out_sim, out_sim_peak);
    nplot = nonneg_plot(out_sim, out_sim_peak);
% 
%     splot = state_plot_3(out, out_sim, out_sim_peak);
% %     
else
    splot = state_plot_N(out, out_sim);
    nplot = nonneg_plot(out_sim);
%     splot = state_plot_2(out, out_sim);
%     splot3 = state_plot_3(out, out_sim);
    
end
