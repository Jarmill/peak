function [fig] = state_plot_2(out, out_sim, out_sim_peak)
%STATE_PLOT_2 Plot states of system trajectories if x has 2 states
%out: information about the recovered solution
%out_sim: function evaluation on random sample trajectories
%out_sim_peak{1}: optional argument, function evaluation if peak is found at
%x0


Tmax = out_sim{1}.t(end);
nx = size(out_sim{1}.x, 2);

box_margin = 0.4*[-1, 1];

assert(nx==2);

fig = figure;
clf
for i = 1:length(out_sim)
    t_curr = out_sim{i}.t;
    x_curr = out_sim{i}.x;
    if i == 1
    
        subplot(1, 3, 1)
        axis square
        hold on
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
        title('Phase Plane of Trajectories')
        
        xlabel('x_1')
        ylabel('x_2')
        
        legend('location', 'northwest') 
        
        subplot(1,3, [2,3])
        hold on        
        plot3(out_sim{i}.t, out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
        pbaspect([2 1 1])         %2.25 1 1
        xlabel('time')
        ylabel('x_1')
        zlabel('x_2')
        if out.optimal
            title(['Peak Value for Trajectories = ', num2str(out.peak_val, 3), ' (optimal),  order = ', num2str(out.order)])
        else
            title(['Peak Value for Trajectories = ', num2str(out.peak_val, 3), ', order = ', num2str(out.order)])
        end
        
        
        legend('location', 'northwest') 
        
    %end
    else
        subplot(1,3, 1)
        
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'off');
        

        subplot(1,3, [2,3])
        
        plot3(out_sim{i}.t, out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'off');
        
    end
end


for i = 1:length(out_sim)
    if i == 1
       subplot(1,3, 1)
        
       scatter(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), 100, 'k', 'DisplayName', 'Initial Points');        

       subplot(1,3, [2,3])
        
       scatter3(0, out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), 100, 'k', 'DisplayName', 'Initial Points');
       
    else
       subplot(1,3, 1)
        
       scatter(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), 100, 'k', 'HandleVisibility', 'off');        

       subplot(1,3, [2,3])
        
       scatter3(0, out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), 100, 'k', 'HandleVisibility', 'off');
        
    
    end
end

%implicit curves
syms y [2 1]
syms t

if out.dynamics.time_indep
    %time independent
    vy = out.func.vval(y) - out.peak_val;
    subplot(1,3,1)
    fimplicit(vy, [xlim+box_margin, ylim + box_margin], ':k', 'DisplayName', 'Invariant Set', 'LineWidth', 3);
    vyt = 1e-8*t + vy ;
else
    vyt = out.func.vval(t, y) - out.peak_val;
end    
cy = out.func.cost(y) + out.peak_val;
cyt = 1e-8*(t + sum(y)) + cy;
fimplicit(cy + 1e-8*sum(y), [xlim+box_margin, ylim + box_margin],  '--r', 'DisplayName', 'Cost Bound', 'LineWidth', 2)

subplot(1,3,[2,3])
xlim([0, Tmax])
fimplicit3(vyt, [xlim, ylim+box_margin, zlim + box_margin], 'EdgeColor', 'None','FaceColor', 'k', 'FaceAlpha', 0.3, ...
            'DisplayName', 'Invariant Set')


fimplicit3(cyt, [xlim, ylim+box_margin, zlim + box_margin], 'EdgeColor', 'None','FaceColor', 'r', 'FaceAlpha', 0.3, ...
            'DisplayName', 'Cost Bound')

view(62, 17)



if out.optimal
    %plot the peak functions too
    subplot(1,3, 1)
        
    plot(out_sim_peak{1}.x(:, 1), out_sim_peak{1}.x(:, 2), 'b', 'HandleVisibility', 'off');

    scatter(out.x0(1), out.x0(2), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);        
    scatter(out.xp(1), out.xp(2), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);        
 
    
    subplot(1,3, [2,3])
        
    plot3(out_sim_peak{1}.t, out_sim_peak{1}.x(:, 1), out_sim_peak{1}.x(:, 2), 'b', 'HandleVisibility', 'off');
    
    
    scatter3(0, out.x0(1), out.x0(2), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);        
    if isfield(out, 'tp')
        scatter3(out.tp, out.xp(1), out.xp(2), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);        
    end
             
end

end

