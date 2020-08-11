%spacecraft attitude control system with 1 degree of freedom 

%https://homepages.laas.fr/henrion/papers/safev.pdf
%find maximum angular velocity with hybrid controller
mset clear
rng(300, 'twister')

%states/variables
mpol('x', 2, 1);


%constants
%these values are high. convert to radians instead of degrees?
I = 27500; %inertia
kp = 2475; %controller gains
kd = 19800;
L = 380; %input saturation


% theta_max = 10*pi/180;
% omega_max = 1*pi/180;
theta_max = 15*pi/180;
omega_max = 3*pi/180;


%theta_max = 5*pi/180;
%omega_max = 0.4*pi/180;
Tmax = 50; %final time

%initial set and feasible space
X0 = [x(1)^2 <= theta_max^2; x(2)^2 <= omega_max^2];
Xsupp = [];


%controller
K = -[kp kd]/L;
u = K*x;

%dynamics with saturation
f_lin = [x(2); u*L/I];  %linear regime
f_up  = [x(2); L/I];        %upper saturation
f_low = [x(2); -L/I];       %lower saturation

X_lin = [u^2 <= 1];
X_up  = [u >= 1];
X_low = [u <= -1];

f = {f_lin, f_up, f_low};
X = {X_lin, X_up, X_low};


%objective to maximize
objective = x(2)^2;        %maximum angular velocity
% objective = u^2;

%formulate problem
p_opt = peak_options;
p_opt.var.x = x;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.state_init = X0;
p_opt.state_supp = Xsupp;
%p_opt.box = [-1, 1; -1, 1; -4, 4];
p_opt.box = [1; 1];
%p_opt.box = [theta_max; omega_max];
p_opt.scale = 0;

p_opt.obj = objective;
%p_opt.Tmax = 20;
order = 4;
out = peak_estimate(p_opt, order);
ang_max = sqrt(-out.peak_val);

%visualize

%Nsample = 100;
Tmax_sim = 30;
Nsample = 50;
sampler = @() (2*rand(2, 1)-1) .* [theta_max; omega_max];

out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim);

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
