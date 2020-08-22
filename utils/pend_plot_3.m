function [fig] = pend_plot_3(out, out_sim, out_sim_peak)
%STATE_PLOT_2 Plot states of system trajectories if x has 2 states
%out: information about the recovered solution
%out_sim: function evaluation on random sample trajectories
%out_sim_peak{1}: optional argument, function evaluation if peak is found at
%x0


Tmax = out_sim{1}.t(end);
nx = size(out_sim{1}.x, 2);

box_margin = 0.4*[-1, 1];

assert(nx==3);

fig = figure;
clf
for i = 1:length(out_sim)
    t_curr = out_sim{i}.t;
    x_curr = out_sim{i}.x;
    if i == 1
    
        
        axis square
        hold on
        plot3(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), out_sim{i}.x(:, 3),  'c', 'DisplayName', 'Trajectories');
        title('Phase Plane of Trajectories')
        
        xlabel('cos(\theta)')
        ylabel('sin(\theta)')
        zlabel('\omega')
        
        legend('location', 'northwest')       
               
    %end
    else
        
        
        plot3(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), out_sim{i}.x(:, 3),  'c', 'HandleVisibility', 'Off');        
    end
end

%initial conditions
for i = 1:length(out_sim)
    if i == 1
       scatter3(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), out_sim{i}.x(1, 3), 100, 'k', 'DisplayName', 'Initial Points');       
    else        
       scatter3(out_sim{i}.x(1, 1), out_sim{i}.x(1, 2), out_sim{i}.x(1, 3), 100, 'k', 'HandleVisibility', 'Off');           
    end
end

scatter3(-1, 0, 0, 300, 'k', 'DisplayName', 'Final Point', 'LineWidth', 2);           

xlim([-1, 1])
ylim([-1, 1])

% %implicit curves
% syms y [3 1]
% syms t
% 
% MD = 60;
% 
% if ~isfield(out, 'tp')
%     %time independent
%     vy = out.func.vval(y) - out.peak_val;    
%     fimplicit3(vy, [xlim+box_margin, ylim + box_margin, zlim + box_margin], 'EdgeColor', 'None','FaceColor', 'k', 'FaceAlpha', 0.3, ...
%             'DisplayName', 'Invariant Set', 'MeshDensity', MD)
%     %vyt = 1e-8*t + vy ;
% else
%     %vyt = out.func.vval(t, y) - out.peak_val;
% end    
% cy = out.func.cost(y) + out.peak_val;
% 
% fimplicit3(cy + 1e-8*sum(y), [xlim+box_margin, ylim + box_margin, zlim + box_margin], 'EdgeColor', 'None','FaceColor', 'r', 'FaceAlpha', 0.3, ...
%             'DisplayName', 'Cost Bound', 'MeshDensity', MD)
% 

view(62, 17)

if nargin == 2
    %plot the peak functions too
    
    plot3(out_sim_peak{1}.x(:, 1), out_sim_peak{1}.x(:, 2), out_sim_peak{1}.x(:, 3), 'b', 'DisplayName', 'Peak Traj.', 'LineWidth', 2);
    scatter3(out_sim_peak{1}.x(1, 1), out_sim_peak{1}.x(1, 2), out_sim_peak{1}.x(1, 3), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);
    
    %initial condition
%     scatter3(out.x0(1), out.x0(2), out.x0(3), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);        
%     scatter3(out.xp(1), out.xp(2), out.xp(3), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);        
end

end

