%% peak impulse response for a state space system
%from the Chesi papers

%starting to turn this into functions

%right now deterministic system
%then move to uncertainty



SOLVE = 1;
PLOT = 1;
MD = 200; %mesh density for 3d
Tmax = Inf;
Tmax_plot = 20;

if SOLVE
    
order = 2;
ranktol = 0.05;

rng(890);

nsys = 1;    
if nsys == 1
    %DC motor example
    %constant Jm
    bm = 0.2;
    Kt = 1;
    La = 0.5;
    Ra = 1;
    Ke = 0.5;
    Jm = 1;
    
    A = [0, 1, 0; 0, -bm/Jm, Kt/Jm; 0, -Ke/La, -Ra/La];
    B = [0; 0; 1/La];
    C = [1, 0, 0];
    
    sys = ss(A, B, C, 0);
    
elseif nsys == 2
    %random system
    n = 6;
    sys = rss(n);
else
    %second order system
    %true upper bound is p* = 0.645
    sys = ss([0 1; -0.5 -1], [0; 1], [1 0], 0);
    %sys = ss([0 1; -0.5 -1], [0; 1], [2 1], 0);
end

n = size(sys.A, 1);

[peak_val, opt] =  peak_impulse_siso(sys.A, sys.B, sys.C, order, ranktol);

end

if PLOT
    figure(1)
clf
xtraj = struct;
[xtraj.y, xtraj.t, xtraj.x] = impulse(sys, Tmax_plot);
margin = 0.2;
xmin = min(xtraj.x) - 0.3;
xmax = max(xtraj.x) + 0.3;
trange = xtraj.t([1, end]);

if (n==2) || (n==3)
subplot(1,2,1)
hold on
if n == 3

    plot3(xtraj.x(:, 1), xtraj.x(:, 2), xtraj.x(:, 3), 'DisplayName', 'state');
    scatter3(sys.B(1), sys.B(2), sys.B(3), 100, 'ok', 'DisplayName', 'Initial Condition', 'LineWidth', 2)
    %fi = fimplicit3(opt.p + opt.obj, [xmin(1), xmax(1), xmin(2), xmax(2),xmin(3), xmax(3)], 'MeshDensity',MD, ...
    %        'EdgeColor', 'None','FaceColor', 'k', 'FaceAlpha', 0.3, ...
    %        'DisplayName', 'Safety Contour');
    if opt.rankp == 1
        scatter3(opt.xp(1), opt.xp(2), opt.xp(3), 100, '*k', 'DisplayName', 'Peak Value', 'LineWidth', 2)
    end
    zlabel('x_3')
else
    %n = 2
    plot(xtraj.x(:, 1), xtraj.x(:, 2), 'DisplayName', 'state', 'LineWidth', 3);
    scatter(sys.B(1), sys.B(2), 100, 'ok', 'DisplayName', 'Initial Condition', 'LineWidth', 2)
    
    p_neg = line_range([sys.C, -peak_val], [xmin(1), xmax(1)], [xmin(2), xmax(2)]);
    p_pos = line_range([sys.C, peak_val], [xmin(1), xmax(1)], [xmin(2), xmax(2)]);
    
    if opt.rankp == 1
        scatter(opt.xp(1), opt.xp(2), 100, '*k', 'DisplayName', 'Peak Value', 'LineWidth', 2)
    end
    if ~isempty(p_neg)
        plot([p_neg{1}(1), p_neg{2}(1)],  [p_neg{1}(2), p_neg{2}(2)], 'r--', 'Linewidth', 3, 'DisplayName', 'Peak Certification')
    end
    if ~isempty(p_pos)
        plot([p_pos{1}(1), p_pos{2}(1)],  [p_pos{1}(2), p_pos{2}(2)], 'r--', 'Linewidth', 3, 'HandleVisibility','off')
    end

    
    fi2 = fimplicit(opt.p  + opt.obj, [xmin(1), xmax(1), xmin(2), xmax(2)], ...
            ':k', 'DisplayName', 'Peak Contour', 'LineWidth', 3, 'MeshDensity', 150);    
end
xlabel('x_1')
ylabel('x_2')
title(['Impulse Response State, order = ', num2str(order)])
axis tight
legend('location', 'northeast')
    

subplot(1,2,2)
end
hold on
title(['Impulse Response Output, peak = ', num2str(peak_val, 3)])
plot(xtraj.t, xtraj.y, 'Linewidth', 3, 'DisplayName', 'Impulse Response')
plot(trange, peak_val*[1,1], 'r--', 'Linewidth', 3, 'DisplayName', 'Peak Certification')
plot(trange, -peak_val*[1,1], 'r--', 'Linewidth', 3, 'HandleVisibility','off')
plot(trange, [0, 0], ':k', 'HandleVisibility','off')
legend('location', 'southeast')

xlabel('t')
ylabel('y')


%% Nonnegative values
figure(2)
clf
%get data along impulse response trajectory 
xcell = num2cell(xtraj.x, 1);
xtraj.v = opt.pval(xcell{:});
xtraj.Lv = opt.Lpval(xcell{:});
xtraj.cost = (xtraj.x * sys.C').^2;
    

   subplot(3,1,1)
    
    hold on        
    plot(xtraj.t, xtraj.v + opt.obj, 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
    
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')

    hold off
    title('Safety Function along Trajectories')
    xlabel('time')
    ylabel('v(x) - \gamma')
    legend('location', 'northwest')
    
    subplot(3,1,2)
    hold on
    plot(xtraj.t, xtraj.Lv, 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')

    hold off
    title('Change in Safety Function along Trajectories')
    xlabel('time')
    ylabel('L_f v(x)')
    legend('location', 'northwest')
    
    subplot(3,1,3)    
    hold on
    
    p_comp = -xtraj.cost - xtraj.v;
    plot(xtraj.t, p_comp, 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
    
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')
    hold off
    
    
    title('Comparision with Objective Function')
    xlabel('time')
    ylabel('-cost(x) - v(x,t)')
    legend('location', 'northwest')
 


end