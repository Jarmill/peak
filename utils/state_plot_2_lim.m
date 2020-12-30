function [fig] = state_plot_2_lim(out, out_sim, out_sim_peak)
%STATE_PLOT_2_LIM Plot states of system trajectories if x has 2 states
% Limited version, does not include (t, x) plot
%out: information about the recovered solution
%out_sim: function evaluation on random sample trajectories
%out_sim_peak{1}: optional argument, function evaluation if peak is found at
%x0


Tmax = out_sim{1}.t(end);
nx = size(out_sim{1}.x, 2);

box_margin = 2;


assert(nx==2);

fig = figure;
clf
for i = 1:length(out_sim)
    t_curr = out_sim{i}.t;
    x_curr = out_sim{i}.x;
    if i == 1
        axis square
        hold on
        if out.dynamics.discrete
            scatter(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), '.c', 'DisplayName', 'Trajectories');
        else
            plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
        end
        title(['Phase Plane of Trajectories, order = ', num2str(out.order)])
        
        xlabel('x_1')
        ylabel('x_2')
        
        legend('location', 'northwest') 
             
    %end
    else
        
        if out.dynamics.discrete
            scatter(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), '.c', 'HandleVisibility', 'off');        
        else
            plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'off');
        end
        
        
        
    end
end

%Initial conditions
for i = 1:length(out_sim)
    if i == 1
        
       scatter(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), 50, 'k', 'DisplayName', 'Initial Points');        

    else
        
       scatter(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), 50, 'k', 'HandleVisibility', 'off');        

    
    end
end

%implicit curves
syms y [2 1]
syms t
MD = 60;
MD = 100;
xlim(stretch(xlim, box_margin))
ylim(stretch(ylim, box_margin))

% PLOT_V = isempty(out.var.w) && isempty(out.var.d) && isempty(out.var.b);
PLOT_V = isempty(out.var.w);

if PLOT_V
    if out.dynamics.time_indep
        %time independent
        vy = -out.func.vval(t, y, []) + out.peak_val;

%         if out.dynamics.Tmax < Inf
%             vy = vy - out.dynamics.Tmax * out.alpha_time;
%         end
        
        fimplicit(vy, [xlim, ylim], ':k', 'DisplayName', 'Invariant Set', 'LineWidth', 3);
    end    
end

nobj = length(out.func.beta);
cyt = cell(nobj, 1);
for i = 1:nobj % multiple objectives
    curr_cost = @(y) out.func.cost_all{i}(y);


    cy = curr_cost(y) - out.peak_val;
    cyt{i} = 1e-8*(t + sum(y)) + cy;
    if i == 1
        fimplicit(cy + 1e-8*sum(y), [xlim, ylim],  '--r', 'DisplayName', 'Cost Bound', 'LineWidth', 2)
    else
        fimplicit(cy + 1e-8*sum(y), [xlim, ylim],  '--r', 'HandleVisibility', 'Off', 'LineWidth', 2)
    end
end


if out.optimal && (nargin == 3)
    %plot the peak functions too
    npeak_traj = length(out_sim_peak);
    for k = 1:npeak_traj
        if k == 1
            if out.dynamics.discrete
                scatter(out_sim_peak{k}.x(:, 1), out_sim_peak{k}.x(:, 2), 60, '.b', 'HandleVisibility', 'off', 'Linewidth', 2);
            else
                plot(out_sim_peak{k}.x(:, 1), out_sim_peak{k}.x(:, 2), 'b', 'HandleVisibility', 'off', 'Linewidth', 2);
            end
            scatter(out.x0(1, k), out.x0(2, k), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);        
            scatter(out.xp(1, k), out.xp(2, k), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);        
           
        else
            
            if out.dynamics.discrete
                scatter(out_sim_peak{k}.x(:, 1), out_sim_peak{k}.x(:, 2), 60, '.b', 'HandleVisibility', 'off');
            else
                plot(out_sim_peak{k}.x(:, 1), out_sim_peak{k}.x(:, 2), 'b', 'HandleVisibility', 'off', 'Linewidth', 2);
            end
            scatter(out.x0(1, k), out.x0(2, k), 200, 'ob', 'HandleVisibility', 'off', 'LineWidth', 2);        
            scatter(out.xp(1, k), out.xp(2, k), 200, '*b', 'HandleVisibility', 'off','LineWidth', 2);        

        end
    end
             
end

end

