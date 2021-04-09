%test out arguments of peak_estimate routine
%clear
SOLVE = 0;
SAMPLE = 0;
PLOT = 0;
V_PLOT = 1;

%% Setup
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
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

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%objective to maximize
objective = -x(2);
%objective = -x(2) - x(1);

% objective = [-x(2); -x(1); -x(1) - x(2)];
% objective = [-x(2); -x(1)];
%
p_opt = peak_options;
p_opt.var.x = x;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;
p_opt.dynamics.discrete = 0;

Tmax_sim = 10;
p_opt.Tmax = Tmax_sim;

p_opt.box = 2.5;
p_opt.scale = 0;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

%% Solving
if SOLVE
    order_max = 5;
    out = cell(order_max, 1);
    for i = 1:order_max
        out{i} = peak_estimate(p_opt, i);
    end
    peak_val = cellfun(@(o) o.peak_val, out);
    optimal = cellfun(@(o) o.optimal, out);

    save('flow_min_x2.mat', 'out', 'peak_val', 'optimal')
end

%% Sampling
if SAMPLE
    load('flow_min_x2.mat')
    Nsample = 60;
    % Nsample = 20;
    s_opt = sampler_options;
    s_opt.sample.x = @() sphere_sample(1, 2)'*R0 + C0;
    s_opt.Tmax = Tmax_sim;
    s_opt.parallel = 0;

    tic
    out_sim = sampler(out{order_max}.dynamics, Nsample, s_opt);
    sample_time = toc;
    if (out{order_max}.optimal == 1)
        s_opt.sample.x = @() out{order_max}.x0;
        out_sim_peak = sampler(out{order_max}.dynamics, 1, s_opt);
    end
    
    save('flow_min_x2.mat', 'out', 'peak_val', 'optimal', 'out_sim', 'out_sim_peak')

end

%% Plotting
if PLOT
    load('flow_min_x2.mat')
%     save('flow_min_x2.mat', 'out', 'peak_val', 'optimal')


    BOUNDS = 0;
    RECOVER = 1;

    figure(1)
    clf
    hold on
    
    theta = linspace(0, 2*pi, 100);
    circ = [cos(theta); sin(theta)];
    %half_theta = linspace(pi/2, 3*pi/2, 100);
    %half_circ = [cos(half_theta); sin(half_theta)];
    
    %initial set    
    X0 = C0 + circ*R0;
    
    plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3)
    
    for i = 1:length(out_sim)
        if i == 1
            plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
        else
            plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
        end
    end
    
    if BOUNDS
        if RECOVER
            
            plot(xlim, [1, 1]*-peak_val(end), ':r', 'LineWidth', 2)
            plot(out_sim_peak{1}.x(:, 1), out_sim_peak{1}.x(:, 2), 'b', 'HandleVisibility', 'off', 'Linewidth', 2);

            scatter(out{end}.x0(1), out{end}.x0(2), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);        
            scatter(out{end}.xp(1), out{end}.xp(2), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);  
        else
            for i = 1:length(peak_val)
                plot(xlim, [1, 1]*-peak_val(i), ':r', 'LineWidth', 2)
            end
        end
    else
        quiver(-0.25, 0.75, 0, -2, 'r', 'LineWidth', 2)
    end
    
    xlim([-0.5, 3])
    ylim([-3, 1.5])
    pbaspect([diff(xlim), diff(ylim), 1])
    axis off 
%     text(
    
end
if V_PLOT
    load('flow_min_x2.mat')
    figure(2)
    clf
    title('Auxiliary Function Comparision', 'FontSize', 16)
    hold on
    for i = 2:length(peak_val)
        v_peak = out{i}.func.vval(out_sim_peak{1}.t', out_sim_peak{1}.x', []);
        plot(out_sim_peak{1}.t, v_peak, 'DisplayName', ['v(t, x) at order ', num2str(i)])
    end    
    
    plot(out_sim_peak{1}.t, out_sim_peak{1}.cost, ':r', 'LineWidth', 2, 'DisplayName', 'Objective Curve')
    scatter(out{end}.tp, -out{end}.xp(2), 200, 'ob', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);        
    legend('location', 'southwest')
    xlabel('time $t$', 'interpreter', 'latex', 'FontSize', 14)
    ylabel('objective $-x_2$', 'interpreter', 'latex','FontSize', 14)
end