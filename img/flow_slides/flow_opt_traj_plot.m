load('flow_min_x2.mat')

[m, i] = max(out_sim_peak{1}.cost);

x_pre = out_sim_peak{1}.x(1:i, :);
x_post= out_sim_peak{1}.x((i+1):end, :);

figure(4)
clf
hold on
% plot(out_sim_peak{1}.x(:, 1), out_sim_peak{1}.x(:, 2), 'b', 'HandleVisibility', 'off', 'Linewidth', 2);
plot(x_pre(:, 1), x_pre(:, 2), 'b', 'Linewidth', 2)
plot(x_post(:, 1), x_post(:, 2), '-.b', 'Linewidth', 2)
plot(xlim, [1, 1]*-peak_val(end), ':r', 'LineWidth', 2)
scatter(out{end}.x0(1), out{end}.x0(2), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);        
scatter(out{end}.xp(1), out{end}.xp(2), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);  

%xlim
pbaspect([diff(xlim), diff(ylim), 1])
% axis off