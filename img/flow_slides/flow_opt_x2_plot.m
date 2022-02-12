% load('flow_min_x2.mat')
load('flow_distance.mat')

% [m, i] = max(out_sim_peak{1}.cost);

% x_pre = out_sim_peak{1}.x(1:i, :);
% x_post= out_sim_peak{1}.x((i+1):end, :);

figure(4)
clf
hold on
for i = 1:length(out_sim)
    if i == 1
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
    else
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
    end
end
plot(out_sim_peak.x(:, 1), out_sim_peak.x(:, 2), 'b', 'HandleVisibility', 'off', 'Linewidth', 2);
% plot(x_pre(:, 1), x_pre(:, 2), 'b', 'Linewidth', 2)
% plot(x_post(:, 1), x_post(:, 2), '-.b', 'Linewidth', 2)
% plot(xlim, [1, 1]*-peak_val(end), ':r', 'LineWidth', 2)
% scatter(out{end}.x0(1), out{end}.x0(2), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);        
% scatter(out{end}.xp(1), out{end}.xp(2), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);  

theta = linspace(0, 2*pi, 250);
X0 = [cos(theta)*0.4+1.5; sin(theta)*0.4];
plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')

%xlim
       xlim([-0.5, 2.5])
    ylim([-1, 1.25])
pbaspect([diff(xlim), diff(ylim), 1])
% axis off