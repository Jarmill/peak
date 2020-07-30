function [fig] = state_plot_N(out, out_sim, out_sim_peak)
%STATE_PLOT_N Plot states of system trajectories if x has >3 states
%out: information about the recovered solution
%out_sim: function evaluation on random sample trajectories
%out_sim_peak{1}: optional argument, function evaluation if peak is found at
%x0

Tmax = out_sim{1}.t(end);
nx = size(out_sim{1}.x, 2);

nplt = nx;

fig = figure;
clf
for i = 1:length(out_sim)
    t_curr = out_sim{i}.t;
    x_curr = out_sim{i}.x;
    for k = 1:(nplt)
        subplot(nplt, 1, k)
        hold on
        
        
        if i == 1
            plot(t_curr, x_curr(:, k), 'c','DisplayName', 'Trajectories')
            if k == 1
                title({['Peak Estimate = ', num2str(out.peak_val), ' , order = ', num2str(out.order)], ['State ', num2str(k), ' vs. Time']})
            else
            title(['State ', num2str(k), ' vs. Time'])
            end
            xlabel('time')
            ylabel(['x', num2str(k)])                
            legend('location', 'east')
            
        else
            plot(t_curr, x_curr(:, k), 'c','HandleVisibility', 'off')                
        end
        if i == length(out_sim)
            plot(xlim, [0, 0], ':k', 'HandleVisibility', 'off')
        end
    end       
end

if out.optimal
    for k = 1:(nplt)
        subplot(nplt, 1, k)
        hold on
        plot(out_sim_peak{1}.t, out_sim_peak{1}.x(:, k), 'b', 'LineWidth', 3, 'DisplayName', 'Peak Traj.')       
        
        if k == 1
            title({['Peak Estimate = ', num2str(out.peak_val), ' (optimal), order = ', num2str(out.order)], ['State ', num2str(k), ' vs. Time']})
        end
    end       
    
end

end

