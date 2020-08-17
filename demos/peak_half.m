% Try the peak estimation on a vector field
% Dynamics from 'On Analysis and Synthesis of Safe Control Laws'
% by Anders Rantzer and Stephen Prajna

%Author: Jared Miller 8/17/20

mset clear
rng(300, 'twister')
SOLVE = 1;
PLOT = 1;

ranktol = 5e-3;

%dynamics
mpol('x', 2, 1);

%support
Xsupp = [];

%dynamics
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
X = [];


%initial set
%C0 = [1.2; 0];
C0 = [1.5; 0];
R0 = 0.4;

%C0 = [0.8; 0];
%R0 = 0.1;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);


%unsafe set (objective to maximize)
%half circle: intersection of disk and half-plane through circle center

%circle
% Cu = [0; -1];
% Cu = [2.5; 0];
Cu = [3; 1];
%Cu = [2; 0]; %unsafe with 7pi/4
%Cu = [2; -0.2];
Ru = 0.4;
cost_circ = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;


%line
%tilt angle
%theta_c = 5*pi/4;
% theta_c = 3*pi/4;
% theta_c = 3*pi/2;    
theta_c = 7*pi/4;


w_c = [cos(theta_c); sin(theta_c)];
cost_line = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

objective = [cost_circ; cost_line];

if SOLVE
    p_opt = peak_options;
    p_opt.var.x = x;

    p_opt.state_supp = Xsupp;
    p_opt.state_init = X0;

    p_opt.dynamics = struct;
    p_opt.dynamics.f = f;
    p_opt.dynamics.X = X;

    Tmax_sim = 10;
    %p_opt.Tmax = Tmax_sim;

    p_opt.box = 4;
    p_opt.scale = 0;
    %p_opt.R = 6;


    p_opt.rank_tol = 4e-3;
    p_opt.obj = objective;

    order = 4;
    out = peak_estimate(p_opt, order);
    peak_val = out.peak_val;
end

if PLOT

    mu = 1;

    Nsample = 100;
%     Nsample = 20;
    sampler = @() circle_sample(1)'*R0 + C0;

    out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, 0, @ode45);

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
    
    theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
    circ_half = [cos(theta_half_range); sin(theta_half_range)];
    Xu = Cu + circ_half* Ru;
    
    subplot(1, 3, 1)
    patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none',...
        'DisplayName', 'Unsafe Set')
    
    
    
    % 
end