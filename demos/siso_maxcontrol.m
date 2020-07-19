SOLVE = 1;
SAMPLE = 1;
PLOT = 1;


% %set up system
%DC motor
% bm = 0.2;
% Kt = 1;
% La = 0.5;
% Ra = 1;
% Ke = 0.5;
% Jm = 1;
% 
% 
% %modifications
% Ra = 0.1;
% Jm = 0.6;
% Kt = 10;
% 
% A = [0, 1, 0; 0, -bm/Jm, Kt/Jm; 0, -Ke/La, -Ra/La];
% B = [0; 0; 1/La];
% C = [1, 0, 0];
% 
% sys = ss(A, B, C, 0);
% 
%weighting
%Q = eye(3);
%R = 15;

%Umich aircraft pitch http://ctms.engin.umich.edu/CTMS/index.php?example=AircraftPitch&section=SystemModeling
A = [-0.313 56.7 0; -0.0139 -0.426 0; 0 56.7 0];
B = [0.232; 0.0203; 0];
C = [0 0 1];
D = [0];
sys = ss(A,B,C,D);


%form controller

%weighting
Q = 2* (C'*C);
R = 1;
%X0 = [0.3; 0; 0.1];
%X0 = [0.5; -0.01; 0.3];
X0 = [0.5; -0.02; 0.7];
%X0 = [0.5; -0.2; 1];

[K, S, e] = lqr(sys, Q, R);

%input to impulse-siso
%X0 = [2; 1; 1];

%X0 = [1; 0; 0];
%X0 = [1; 3; 1];
%X0 = [0; 0; 4];
if SOLVE
    Ac = A - B * K;
    Bc = X0;
    Cc = K;

    order = 2;

    [peak_val, opt] = peak_impulse_siso(Ac, Bc, Cc, order);
end

%% Simulations
if SAMPLE
    Tmax = 20;  %Maximum time to plot
  
    xtraj = struct;
    [xtraj.t, xtraj.x] = ode45(@(t, x) Ac*x, [0, Tmax], Bc);  
    xtraj.u = -xtraj.x * K';
    xcell = num2cell(xtraj.x, 1);
        
    xtraj.v = opt.pval(xcell{:});   
    xtraj.Lv(:) = opt.Lpval(xcell{:});
    %xtraj.cost = (xtraj.x * Cc').^2;
    xtraj.cost = xtraj.u.^2;

    
    %evaluate functions along trajectories
end


if PLOT
    coord_lim = [-1, 1];
    Nsample = 1;
    nsys = 1;
    
    FS = 14;
    
    %% States and Contours
    figure(1)
    clf
    subplot(3,1,1)
    hold on
    plot(xtraj.t, xtraj.x(:, 1), 'DisplayName', 'State 1')
    plot(xtraj.t, xtraj.x(:, 2), 'DisplayName', 'State 2')
    plot(xtraj.t, xtraj.x(:, 3), 'DisplayName', 'State 3')
    plot(xtraj.t, xtraj.u, 'k', 'LineWidth', 2, 'DisplayName', 'Control')
    
    plot([0, Tmax], [peak_val, peak_val], 'r--', 'LineWidth', 2, 'DisplayName', 'Control Bound')
    plot([0, Tmax], -[peak_val, peak_val], 'r--', 'LineWidth', 2, 'HandleVisibility', 'Off')
    
    legend('location', 'northeast')
    title('States and Control Effort')
    hold off
    
    
    subplot(3, 1, [2,3])
    
    hold on
    for i = 1:Nsample
        if i == 1
            plot3(xtraj.x(:, 1), xtraj.x(:, 2), xtraj.x(:, 3), 'color', 0.6*[1,1,1], 'DisplayName', 'Trajectories')
        else
            plot3(xtraj.x(:, 1), xtraj.x(:, 2), xtraj.x(:, 3), 'color', 0.6*[1,1,1], 'HandleVisibility', 'Off')
        end
    end

    title(['Linear System, order = ', num2str(order), ', peak estimate = ', num2str(peak_val, 4)], 'FontSize', FS)
    
    xlim(coord_lim)
    ylim(coord_lim/30)
    zlim(coord_lim)
    xlabel('x_1')
    ylabel('x_2')
    zlabel('x_3')
    scatter3(Bc(1), Bc(2),Bc(3), 200, 'ok', 'DisplayName', 'Initial Condition', 'LineWidth', 2)
    %p_neg = line_range([C, -opt.peak_val], [-2, 2], [-2, 2]);
    %p_pos = line_range([C, opt.peak_val],  [-2, 2], [-2, 2]);
    if opt.rankp == 1
        scatter3(opt.xp(1), opt.xp(2), opt.xp(3), 100, '*k', 'DisplayName', 'Peak Value', 'LineWidth', 2)
    end
    
    
    if K(3)==0
        %The power of patches! (projective geometry, assuming that x != 0)
        %make this actually robust
        zplane = [1 -1 -1 1]; % Generate data for x vertices
        yplane = [1 1 -1 -1]; % Generate data for y vertices
        xplane_pos = -1/K(1)*(K(3)*zplane + K(2)*yplane + peak_val); % Solve for z vertices data
        xplane_neg = -1/K(1)*(K(3)*zplane + K(2)*yplane - peak_val); % Solve for z vertices data
        
    else
        xplane = [1 -1 -1 1]; % Generate data for x vertices
        yplane = [1 1 -1 -1]; % Generate data for y vertices
        zplane_pos = -1/K(3)*(K(1)*xplane + K(2)*yplane + peak_val); % Solve for z vertices data
        zplane_neg = -1/K(3)*(K(1)*xplane + K(2)*yplane - peak_val); % Solve for z vertices data
        
    end
    
    patch('XData',xplane, 'YData',yplane, 'ZData', zplane_pos, ...
        'FaceColor', 'r', 'FaceAlpha', 0.3, 'DisplayName', 'Maximum Control')
    patch('XData',xplane, 'YData',yplane, 'ZData', zplane_neg, ...
        'FaceColor', 'r', 'FaceAlpha', 0.3, 'HandleVisibility', 'Off')
    
%     if ISOSURFACE
%         MD = 210;
%         fi = fimplicit3(opt.p + opt.obj, [coord_lim, coord_lim/30, coord_lim], 'MeshDensity',MD, ...
%            'EdgeColor', 'None','FaceColor', 'k', 'FaceAlpha', 0.3, ...
%            'DisplayName', 'Safety Contour');
%     end
    %if ~isempty(p_neg)
    %    plot([p_neg{1}(1), p_neg{2}(1)],  [p_neg{1}(2), p_neg{2}(2)], 'r--', 'Linewidth', 3, 'DisplayName', 'Peak Certification')
    %end
    %if ~isempty(p_pos)
    %    plot([p_pos{1}(1), p_pos{2}(1)],  [p_pos{1}(2), p_pos{2}(2)], 'r--', 'Linewidth', 3, 'HandleVisibility','off')
    %end

    %fi2 = fimplicit(opt.p  + opt.obj, [-2, 2, -2, 2], ...
    %            ':k', 'DisplayName', 'Peak Contour', 'LineWidth', 3, 'MeshDensity', 150);    
    hold off
    axis square
    legend('location', 'southwest')
    
    %% Values
    figure(2)
    clf
    tiledlayout(1+2, 1);
    
%     subplot(3,1,1)
    nexttile
    hold on       
    for i = 1:Nsample
        if i == 1
            plot(xtraj.t, xtraj.v  + opt.obj, 'c','DisplayName', 'Trajectories')
        else
            plot(xtraj.t, xtraj.v + opt.obj, 'c', 'HandleVisibility','off')
        end  
    end
    
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')

    hold off
    title('Safety Function along Trajectories', 'FontSize', FS)
    xlabel('time')
    ylabel('v(x) - \gamma')

    
%     subplot(3,1,2)
    
    for j = 1:nsys
        nexttile
        hold on
        for i = 1:Nsample
            if i == 1
                plot(xtraj.t, xtraj.Lv, 'c','DisplayName', 'Trajectories')
            else
                plot(xtraj.t, xtraj.Lv, 'c', 'HandleVisibility','off')
            end  
        end

        %plots

        plot(xlim, [0, 0], ':k',  'HandleVisibility','off')

        hold off
        title(['Change in Safety Function for System ', num2str(j)], 'FontSize', FS)
        xlabel('time')
        ylabel(['L_{f', num2str(j), '}v(x)'])
    
    end
%     subplot(3,1,3)  
    nexttile
    hold on
    
    for i = 1:Nsample
        p_comp = -xtraj.cost - xtraj.v;
        
        if i == 1
            plot(xtraj.t, p_comp, 'c','DisplayName', 'Trajectories')
        else
            plot(xtraj.t, p_comp, 'c', 'HandleVisibility','off')
        end  
    end

    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')
    hold off
    
    
    title('Comparision with Objective Function', 'FontSize', FS)
    xlabel('time')
    ylabel('-cost(x) - v(x,t)')

    
end