function [fig] = nonneg_plot(out_sim, out_sim_peak)
%NONNEG_PLOT Plot nonnegative functions (value, Lie derivatives, cost
%comparision)
%out_sim: function evaluation on random sample trajectories
%out_sim_peak: optional argument, function evaluation if peak is found at
%x0


%if iscell(f)
%    nsys = length(f);
%else
%    nsys = 1;
%end

nplt = size(out_sim{1}.nonneg, 1);

fig = figure;
clf
for i = 1:length(out_sim)
    t_curr = out_sim{i}.t;
    nonneg_curr = out_sim{i}.nonneg;
    for k = 1:(nplt)
        subplot(nplt, 1, k)
        hold on
        
        
        if i == 1
            if k == 1
                title('Value Function')
                plot(t_curr, nonneg_curr(k, :), 'c','DisplayName', 'Trajectories')
                title('Safety Function along Trajectories')
                xlabel('time')
                ylabel('v(x,t) - \gamma')                
                legend('location', 'east')
            elseif k == nplt
                title('Cost Comparision')
                xlabel('time')
                ylabel('-cost(x) - v(x,t)')
                legend('location', 'east')
                plot(t_curr, nonneg_curr(k, :), 'c', 'DisplayName', 'Trajectories')
            else
                title(['Change in Value Function along System ', num2str(k-1)])
                xlabel('time')
                ylabel(['L_{f', num2str(k-1), '} v(x,t)'])
                legend('location', 'east')
                plot(t_curr, nonneg_curr(k, :), 'c', 'DisplayName', 'Trajectories')
            end
            
        else
            plot(t_curr, nonneg_curr(k, :), 'c','HandleVisibility', 'off')                
        end
        if i == length(out_sim)
            plot(xlim, [0, 0], ':k', 'HandleVisibility', 'off')
        end
    end       
end

if nargin == 2  
    for k = 1:(nplt)
        subplot(nplt, 1, k)
        hold on
        plot(out_sim_peak{1}.t, out_sim_peak{1}.nonneg(k, :), 'b', 'LineWidth', 3, 'DisplayName', 'Peak Traj.')       
    end       
end

end

