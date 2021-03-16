%spacecraft attitude control system with 1 degree of freedom 

%theta_max = 20.6899 degrees
%order = 5

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
% objective = x'*x;
objective = x(1)^2;
%objective = x(2)^2;        %maximum angular velocity
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
order = 5;
out = peak_estimate(p_opt, order);
ang_max = sqrt(-out.peak_val);

%visualize


%deal with symmetry
M0_mix = out.M0(2:3, 2:3);
Mp_mix = out.Mp(2:3, 2:3);

rank_sym_0 = rank(M0_mix, p_opt.rank_tol);
rank_sym_p = rank(Mp_mix, p_opt.rank_tol);


%% symmetry
if rank_sym_0==1 && rank_sym_p==1
    %the solution is optimal with symmetry
    %the symmetry structure of the solution is x <=> -x
    %equivariance: f(-x) = -f(x)
    

    
    %find orbit of best initial point x0
    sym_sign0 = sign(M0_mix(2,1));
    x0_abs = sqrt([M0_mix(1,1); M0_mix(2,2)]);
    x0 = [1; sym_sign0].* x0_abs;
%     x0_neg = -x0_pos;
    
    %find orbit of peak achieved point xp
    sym_signp = sign(Mp_mix(2,1));    
    xp_abs = sqrt([Mp_mix(1,1); Mp_mix(2,2)]);
    xp = [1; sym_signp].* xp_abs;
%     xp_neg = -xp_pos;

    %package into the output
    if ~out.dynamics.time_indep
        out.tp = [1 1] * out.tp;
    end
    out.x0 = [x0 -x0];
    out.xp = [xp -xp];
%     out.w = [1 1]*out.w;
    out.optimal = 1;
end
    

%% visualize
%Nsample = 100;
Tmax_sim = 30;
Nsample = 50;
% sampler = @() (2*rand(2, 1)-1) .* [theta_max; omega_max];

s_opt = sampler_options;
% s_opt.sample.x = @() circle_sample(1)'*R0 + C0;
s_opt.sample.x = @() (2*rand(2, 1)-1) .* [theta_max; omega_max];
s_opt.Tmax = Tmax_sim;


% out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim);
tic
out_sim = sampler(out.dynamics, Nsample, s_opt);
sample_time = toc;
if (out.optimal == 1)
%     out_sim_peak = switch_sampler(out.dynamics, out.x0, 2, Tmax_sim);
    s_opt.sample.x = out.x0;
    out_sim_peak = sampler(out.dynamics, 2, s_opt);
    
    nplot = nonneg_plot(out, out_sim, out_sim_peak);

    splot = state_plot_2(out, out_sim, out_sim_peak);
    %     splot = state_plot_N(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
    splot = state_plot_2(out, out_sim);
%     splot = state_plot_N(out, out_sim);
end


%% saturation plot
subplot(1,3,1)
p = line_range([K 1], xlim, ylim);
plot([p{1}(1), p{2}(1)], [p{1}(2), p{2}(2)], ':k', 'DisplayName', 'Saturation');
p = line_range([K -1], xlim, ylim);
plot([p{1}(1), p{2}(1)], [p{1}(2), p{2}(2)], ':k', 'HandleVisibility', 'Off');
