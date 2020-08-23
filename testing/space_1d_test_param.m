%spacecraft attitude control system with 1 degree of freedom 
%https://homepages.laas.fr/henrion/papers/safev.pdf

%find maximum angular velocity with hybrid controller
mset clear

%order 5
%p = 51.85 degrees


%states/variables
mpol('x', 2, 1);
mpol('w', 1, 1);


%constants
%these values are high. convert to radians instead of degrees?
I = 27500; %inertia
kp = 2475; %controller gains
kd = 19800;
L = 380; %input saturation
Wmax = 0.5;
% Wmax = 0.1;


% theta_max = 10*pi/180;
% omega_max = 1*pi/180;
theta_max = 20*pi/180;
omega_max = 5*pi/180;


%theta_max = 5*pi/180;
%omega_max = 0.4*pi/180;
Tmax = 50; %final time

%initial set and feasible space
X0 = [x(1)^2 <= theta_max^2; x(2)^2 <= omega_max^2];
% X0 = [x(1) <= theta_max; x(2)^2 <= omega_max^2];

Xsupp = [];
Wsupp = [w^2 <= Wmax^2];
%Wsupp = (w == 0);
% Wsupp = [w >= -sqrt(Wmax); w <= 0];

%controller
K = -[kp kd]/L;
u = K*x;

%dynamics with saturation
f_lin = [x(2); u*L/I * (1+w)];  %linear regime
f_up  = [x(2); L/I * (1+w)];        %upper saturation
f_low = [x(2); -L/I * (1+w)];       %lower saturation

X_lin = [u^2 <= 1];
X_up  = [u >= 1];
X_low = [u <= -1];

f = {f_lin, f_up, f_low};
X = {X_lin, X_up, X_low};

% f = f_lin;
% X = [];

%objective to maximize
%objective = x(1)^2;        %maximum angle
% objective = x(2)^2;        %maximum angular velocity
% objective = u^2;
% objective = x'*x;
% objective = (x(1) + 10*x(2));
objective = x(1)^2;

%formulate problem
p_opt = peak_options;
p_opt.var.x = x;
p_opt.var.w = w;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

p_opt.state_init = X0;
p_opt.state_supp = Xsupp;
p_opt.param = Wsupp;
%p_opt.box = [-1, 1; -1, 1; -4, 4];
p_opt.box = 2;
%p_opt.box = [theta_max; omega_max];

p_opt.obj = objective;
%p_opt.Tmax = 20;
order = 5;
out = peak_estimate(p_opt, order);
ang_max = sqrt(-out.peak_val);


%% Deal with symmetry


%deal with symmetry
M0_mix = out.M0(2:3, 2:3);
Mp_mix = out.Mp(2:3, 2:3);

rank_sym_0 = rank(M0_mix, p_opt.rank_tol);
rank_sym_p = rank(Mp_mix, p_opt.rank_tol);

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
    out.w = [1 1]*out.w;
    out.optimal = 1;
end
    

%% visualize

Nsample = 80;
 Tmax_sim = 60;
% Tmax_sim = 30;
% Tmax_sim = 10;
% Nsample = 50;
sampler = @() (2*rand(3, 1)-1) .* [theta_max; omega_max; Wmax];
% 
out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, 10, length(w));
% 
if (out.optimal == 1)
    out_sim_peak = switch_sampler(out.dynamics, [out.x0; out.w], 2, Tmax_sim, 10, length(w));
    nplot = nonneg_plot(out, out_sim, out_sim_peak);

    splot = state_plot_2(out, out_sim, out_sim_peak);
    %     splot = state_plot_N(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
    splot = state_plot_2(out, out_sim);
%     splot = state_plot_N(out, out_sim);
end

subplot(1,3,1)
p = line_range([K 1], xlim, ylim);
plot([p{1}(1), p{2}(1)], [p{1}(2), p{2}(2)], ':k', 'DisplayName', 'Saturation');
p = line_range([K -1], xlim, ylim);
plot([p{1}(1), p{2}(1)], [p{1}(2), p{2}(2)], ':k', 'HandleVisibility', 'Off');
