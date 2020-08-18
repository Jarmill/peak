function [fig] = cost_plot(out, out_sim, out_sim_peak)
%COST_PLOT Plot the objective values along trajectories
%out_sim: function evaluation on random sample trajectories
%out_sim_peak: optional argument, function evaluation if peak is found at x0

fig = figure;
clf
for i = 1:length(out_sim)
    t_curr = out_sim{i}.t;
    cost_curr = out_sim{i}.cost;
    hold on
    
    if i == 1
        plot(t_curr, cost_curr, 'c','DisplayName', 'Trajectories')
        title('Cost along Trajectories')
        xlabel('time')
        ylabel('cost(x)')                
        legend('location', 'southeast')
    else
        plot(t_curr, cost_curr, 'c','HandleVisibility', 'off')                
    end    
    
end

if nargin == 3  
    hold on
    plot(out_sim_peak{1}.t, out_sim_peak{1}.cost, 'b', 'LineWidth', 3, 'DisplayName', 'Peak Traj.')           
    
    if isfield(out, 'tp')
        scatter(out.tp(1), out.peak_val, 300, '*b', 'Linewidth', 2, 'DisplayName', 'Peak Achieved')
    end
end

plot(xlim, [1,1]*out.peak_val, '--r', 'Linewidth', 3, 'DisplayName', 'Cost Bound')


end

