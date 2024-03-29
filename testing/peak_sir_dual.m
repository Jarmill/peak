% Peak Estimation of SIR system

%Author: Jared Miller 6/25/20
SOLVE = 1;
PARAM = 2;

PLOT = 1;

beta0 = 0.4;
% gamma0 = 0.04;
gamma0 = 0.1;
I_max = 0.1;

beta_tol = 0.2;  % 20% uncertainty in infection rate
gamma_tol = 0.1; % 10% uncertainty in removal rate

if SOLVE    
    %mset clear
    %mpol('x', 2, 1);
    x = sdpvar(2, 1);
    if PARAM==2
        %beta and gamma
        w = sdpvar(2, 1);
        Wsupp = struct('ineq',  [beta_tol; gamma_tol].^2 - w.^2);
        beta = beta0*(1 + w(1));
        gamma = gamma0*(1 + w(2));
    elseif PARAM == 1
        %beta only
        w = sdpvar(1, 1);
        Wsupp = struct('ineq', w^2 <= beta_tol^2);
        beta = beta0*(1 + w);
        gamma = gamma0;
    else
        w = [];
        Wsupp = [];
        beta = beta0;
        gamma = gamma0;
    end
    
    %Xsupp = [sum(x) <= 1; x >= 0];
    Xsupp = struct('ineq', [1-sum(x); x]);
    
    %X0 = (x(2) <= I_max);
    X0 = struct('ineq', (I_max - x(2)));
    
    f = [-beta*x(1)*x(2); beta*x(1)*x(2) - gamma*x(2)];
    X = [];
    
    %max number of infected people at any one time
    objective = x(2);
    
    %set up variables
    p_opt = peak_options;
    p_opt.var.x = x;
    p_opt.var.w = w;    

    %dynamics
    p_opt.dynamics = struct;
    p_opt.dynamics.f = f;
    p_opt.dynamics.X = X;
    
    %support sets
    p_opt.state_init = X0;
    p_opt.state_supp = Xsupp;
    p_opt.param = Wsupp;
    
    p_opt.obj = objective;
    
    order = 4;
    out = peak_estimate_dual(p_opt, order);
end

if PLOT       
    rng(300, 'twister')
    %sample from X0 
    %Nsample = 120;
    Nsample = 100;
    Tmax_sim = 30;
    
    rng(33, 'twister')
    if PARAM == 2
        sampler = @() [trap_sample(1, I_max); (2*rand(2, 1)-1).*[beta_tol; gamma_tol]];
    elseif PARAM == 1
        sampler = @() [trap_sample(1, I_max); (2*rand()-1)*beta_tol];
    else
        sampler = @() trap_sample(1, I_max);
    end
    %Xsample = trap_sample(Nsample, I_max);    

    %sample the trajectories
    out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, 10, length(w));
    X0_pts = [0, 1, 1-I_max, 0, 0;
         0, 0, I_max, I_max, 0];
    
    % display the results
    if (out.optimal == 1)        
        out_sim_peak = switch_sampler(out.dynamics, [out.x0; out.w], 1, Tmax_sim, 10, length(w));
        nplot = nonneg_plot(out_sim, out_sim_peak);

        splot = state_plot_2(out, out_sim, out_sim_peak);
        
        %     splot = state_plot_N(out, out_sim, out_sim_peak);
    else
        nplot = nonneg_plot(out_sim);
        splot = state_plot_2(out, out_sim);
    %     splot = state_plot_N(out, out_sim);
    end

    %modify state plot 
    subplot(1, 3, 1)
    plot(X0_pts(1, :), X0_pts(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')
    xlabel('Susceptible')
    ylabel('Infected')
    xlim([0, 1])
    ylim([0, 1.1])
    
    subplot(1, 3, [2,3])
    plot3(zeros(5, 1), X0_pts(1, :), X0_pts(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')
    xlabel('Time')
    ylabel('Susceptible')
    zlabel('Infected')
    ylim([0, 1.1])
    zlim([0, 1.1])
        
end
