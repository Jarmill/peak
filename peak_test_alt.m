% Peak estimation testbed
% Try to find the maximum value of a function p(x)
% on trajectories xdot = f(t,x) starting from a set X0

% Try the peak estimation on a vector field
% Dynamics from 'On Analysis and Synthesis of Safe Control Laws'
% by Anders Rantzer and Stephen Prajna

%Author: Jared Miller 6/22/20
SOLVE = 0;
PLOT = 1;
ISOSURFACE = 0;
MD = 80; %mesh density of isosurface plot


%Enable scaling of all variables through mpol/scale (SCALE = 1)
%or manually scale the time coordinate only (SCALE = 0)
SCALE = 0;

T = 8;
%T = 20;   %final time

%initial set
C0 = [1.5; 0];
R0 = 0.4;

%C0 = [1.5; 1];
%R0 = 0.25;

%unsafe set
Cu = [-1; -1];
Ru = 0.4;

%plotting
%m_low = -1;
m_low = -3;
m_high = 3;

%subspace angle
%theta = 3*pi/2; %(equivalent to maximizing -x(2))
%theta = pi; %max -x(1)
theta = 11*pi/8;
LINE_COST = 1;

%dynamics

%prajna and rantzer
fv = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

%linear decay
%fv = @(t, x) [-x(1); -x(1) - x(2)];

%linear oscillation decay
%fv = @(t, x) [-x(2); x(1) - 0.2*x(2)];

%linear oscillation
%fv = @(t, x) 2*[-x(2); x(1)];

%rusty spring
%fv =  @(t, x) [x(2); -1*x(1) - 0.25*x(1)^3];

%van der pol
%fv = @(t, x) [x(2); (1-x(1)^2)*x(2) - x(1)];

%predator prey
%fv = @(t, x) [(2 * x(1) - 4*x(1)*x(2))/3; x(1)*x(2) - x(2)];

if SOLVE
    %% Parameters for solver    
    mset clear; warning('off','YALMIP:strict')
    mset('yalmip',true);
    mset(sdpsettings('solver', 'mosek'));

    %order = input('order of relaxation ='); d = 2*order;
    order = 4;
    d = 2*order;
        
    R = 5;    %radius to contain dynamics
    %dynamics are time-independent

    %% Formulate Measures and System Properties
    
    %variables and measures
    mpol('t0', 1); mpol('x0', 2); mu0 = meas([t0; x0]); %initial measure
    mpol('t', 1);  mpol('x', 2);  mu  = meas([t; x]);   %occupation measure    
    mpol('tp', 1); mpol('xp', 2); mup = meas([tp; xp]); %peak measure

    %test functions (monomials)
    
    v0 = mmon([t0; x0], d);
    v  = mmon([t; x], d);
    vp = mmon([tp; xp], d);

    %unknown moments of initial measure
    y0 = mom(v0);
    yp = mom(vp);

    %dynamics and cost   
    if SCALE
        f = fv(t, x);
    else
        f = T*fv(t, x);
    end
    
    %Liouville Equation 
    Ay = mom(diff(v, t) + diff(v, x)*f); 
    Liou = Ay + (y0 - yp);

    mom_con = [Liou == 0; mass(mu0)==1];


    %Support Constraints
    X  = (x'*x <= R^2);
    Xp = (xp'*xp <= R^2);
    X0 = ((x0(1)-C0(1))^2 + (x0(2)-C0(2))^2 <= R0^2);
    
    if SCALE
        Tsupp = [t*(T-t) >= 0, tp*(T - tp) >= 0, t0 == 0];
        scale(t, T); scale(tp, T); scale(x, R); scale(x0, R); scale(xp, R);
    else
        Tsupp = [t*(1-t) >= 0, tp*(1 - tp) >= 0, t0 == 0];
    end
        
    supp_con = [X0, X, Xp, Tsupp];
    
    
    %function to optimize
    if LINE_COST                
        costv = @(x) x * [cos(theta); sin(theta)];
    else
        costv = Ru^2 - (x(:, 2)-Cu(1)).^2 - (x(:, 2)-Cu(2)).^2;
    end
    
    cost = costv(xp);
    
    objective = max(cost);

    %% Solve LMI and recover solution
    
    %Input LMI moment problem
    P = msdp(objective, ...
        mom_con, supp_con);

    %solve LMI moment problem    
    [status,obj,m,dual_rec]= msol(P);
    %[model.A, model.b, model.c, model.K, model.b0, model.s] = msedumi(P);
    
    %extract solutions
    gamma_val = -obj;
    radius_out = sqrt(gamma_val + Ru^2);
    M0 = double(mmat(mu0));
    Mp = double(mmat(mup));
    
    %1, t, x1, x2
    M0_1 = M0(1:4, 1:4);
    Mp_1 = Mp(1:4, 1:4);
    
    rank0 = rank(M0_1, 1e-4);
    rankp = rank(Mp_1, 1e-4);
    
    %time moments
    if SCALE
        t0_out = double(mom(t0));        
        tp_out = double(mom(tp));
    else
        t0_out = T*double(mom(t0));        
        tp_out = T*double(mom(tp));
    end
    
    
    %space moments
    x0_out = double(mom(x0));
        
    xp_out = double(mom(xp));    
end

if PLOT
    %% Set and Trajectory Properties
    Nsample = 40;

    %initial and unsafe sets
    Ntheta = 100;
    
    theta_circ = linspace(0, 2*pi, Ntheta);
    circ = [cos(theta_circ); sin(theta_circ)];
    
    %initial set        
    X0_pts = C0 + circ*R0;
    
    %sample from X0 
    rng(33, 'twister')

    Xsample = C0 + circle_sample(Nsample)' * R0;
    
    
    
    %unsafe set    
    Xu = Cu + circ*Ru;            
    Xmargin = Cu + circ*radius_out;
    

    
    figure(1)
    clf
    %% Plot Phase Plane
    subplot(1,3,1)
    hold on
    
    %initial set
    plot(X0_pts(1, :), X0_pts(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')
    
    
    %safety margins and unsafe region
    if LINE_COST
        
        box_pts = line_range([cos(theta) sin(theta) gamma_val], [m_low, m_high], [m_low, m_high]);
        plot([box_pts{1}(1), box_pts{2}(1)], [box_pts{1}(2), box_pts{2}(2)], 'r--', 'Linewidth', 3, 'DisplayName', 'Safety Margin')
        %plot([m_low, m_high], gamma_val*[1, 1], 'r--', 'Linewidth', 3, 'DisplayName', 'Safety Margin')
        
    else
        patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
        plot(Xmargin(1, :), Xmargin(2, :), 'r--', 'Linewidth', 3, 'DisplayName', 'Initial Set')
    end
        
    %trajectories
    xtraj = cell(Nsample, 1);
    for i = 1:Nsample
        xtraj{i} = struct;
        [xtraj{i}.t, xtraj{i}.x] = ode45(fv, [0, T], Xsample(:, i));  
        
        if i == 1
            plot(xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories')
        else
            plot(xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c', 'HandleVisibility','off')
        end
    end 
    
    
    if rank0 == 1 && rankp == 1
        %global optimum was found
        
        %peak trajectory
        xtraj_opt = struct;
        [xtraj_opt.t, xtraj_opt.x] = ode45(fv, [0, T], x0_out);  
        
        plot(xtraj_opt.x(:, 1), xtraj_opt.x(:, 2), 'b', 'LineWidth', 2, 'DisplayName', 'Peak Trajectory')
        
        
        MS = 200;
        
        %initial and optimal points on peak trajectory
        scatter(xp_out(1), xp_out(2), MS, '*b', 'HandleVisibility','off', 'LineWidth', 2)
        scatter(x0_out(1), x0_out(2), MS, 'ob', 'HandleVisibility','off', 'LineWidth', 2)
        
        legend('location', 'northwest')    
    end    
    
    title(['Safety Margin for Trajectories = ', num2str(-gamma_val, 3), ' order = ', num2str(order)])

    xlim([m_low, m_high])
    ylim([m_low, m_high])
    
    hold off
    axis square
    xlabel('x')
    ylabel('y')
    
    %% Plot trajectories in time
    
    subplot(1,3,[2,3])    
    hold on
    
    %initial set
    plot3(zeros(Ntheta, 1), X0_pts(1, :), X0_pts(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')
    
    %lower bound to maximum (safety margin)
    if LINE_COST                
        %plot([box_pts{1}(1), box_pts{2}(1)], [box_pts{1}(2), box_pts{2}(2)], 'Linewidth', 3, 'DisplayName', 'Safety Margin')
        YData = [box_pts{1}(1), box_pts{2}(1), box_pts{2}(1), box_pts{1}(1)];
        ZData = [box_pts{1}(2), box_pts{2}(2), box_pts{2}(2), box_pts{1}(2)];
        
        patch('XData',[0, 0, T, T], 'YData',YData, 'ZData', ZData, ...
            'FaceColor', 'r', 'FaceAlpha', 0.3, 'DisplayName', 'Safety Margin')
    
        %patch('XData',[0, 0, T, T], 'YData',[m_low, m_high, m_high, m_low], ...
        %    'ZData', gamma_val*[1,1,1,1], 'FaceColor', 'r', 'FaceAlpha', 0.3, 'DisplayName', 'Safety Margin')
    else
        [xcyl, ycyl, zcyl] = cylinder(1, 50);
        surf(T*zcyl, Cu(1) + xcyl*Ru, Cu(2) + ycyl*Ru, 'FaceColor', 'r', 'EdgeColor', 'None', 'DisplayName', 'Unsafe Set')
        surf(T*zcyl, Cu(1) + xcyl*radius_out, Cu(2) + ycyl*radius_out, ...
            'FaceColor', 'r', 'EdgeColor', 'None', 'Facealpha', 0.3, 'DisplayName', 'Safety Margin')
    end
    
    %contour in SOS program (safety contour)

    syms tc xc yc;
    %scaling with t
    %vv = monolist([tc/T; xc; yc], d);
    %vv = conj(monolistYalToGlop([tc/T; xc; yc], d));
    if SCALE
        vv = conj(monolistYalToGlop([tc; xc; yc], d));
    else
        vv = conj(monolistYalToGlop([tc/T; xc; yc], d));
    end

    %recovered dual variables from msdp, correspond to the free
    %variables in the sedumi problem
    p = dual_rec'*vv;
    pval = matlabFunction(p);
    
    Lp = diff(p, tc) + [diff(p, xc) diff(p, yc)]*fv(tc, [xc, yc]);
    Lpval = matlabFunction(Lp);

    if ISOSURFACE
        fi = fimplicit3(p  + obj, [0, T, m_low, m_high, m_low, m_high], 'MeshDensity',MD, ...
            'EdgeColor', 'None','FaceColor', 'k', 'FaceAlpha', 0.3, ...
            'DisplayName', 'Safety Contour');

    end
    
    %sample some trajectories in X0
    for i = 1:Nsample        
        if i == 1
            plot3(xtraj{i}.t, xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c','DisplayName', 'Trajectories')
        else
            plot3(xtraj{i}.t, xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c', 'HandleVisibility','off')
        end        
    end  
    
    
    if rank0 == 1 && rankp == 1
        %global optimum was found
        
        %plot peak trajectory
        plot3(xtraj_opt.t, xtraj_opt.x(:, 1), xtraj_opt.x(:, 2), 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')
        
        
        MS = 200;
        
        %initial and optimal points on peak trajectory
        scatter3(tp_out, xp_out(1), xp_out(2), MS, '*b', 'HandleVisibility','off', 'LineWidth', 2)
        scatter3(0, x0_out(1), x0_out(2), MS, 'ob', 'HandleVisibility','off', 'LineWidth', 2)
        
        
        title(['Safety Margin for Trajectories = ', num2str(-gamma_val, 2), ' at time t=', num2str(tp_out, 3)])
    else
        title(['Safety Margin for Trajectories = ', num2str(-gamma_val, 3)])
    end     
    
    legend('location', 'northwest')
    
    xlim([0, T])
    ylim([m_low, m_high])
    zlim([m_low, m_high])
    xlabel('time')
    ylabel('x_1')
    zlabel('x_2')
    pbaspect([2.25 1 1])
    
    
        
    %% compute value functions along trajectories           
    for i = 1:Nsample        
        xtraj{i}.v = pval(xtraj{i}.t, xtraj{i}.x(:, 1),xtraj{i}.x(:, 2));            
        xtraj{i}.x_valid = (sum(xtraj{i}.x.^2, 2) <= R^2);
        xtraj{i}.Lv = Lpval(xtraj{i}.t, xtraj{i}.x(:, 1),xtraj{i}.x(:, 2));                               
        xtraj{i}.cost = costv(xtraj{i}.x);
    end

    if rank0 == 1 && rankp == 1
        %global optimum trajectory
        xtraj_opt.v = pval(xtraj_opt.t, xtraj_opt.x(:, 1), xtraj_opt.x(:, 2));
        xtraj_opt.x_valid = (sum(xtraj_opt.x.^2, 2) <= R^2);
        xtraj_opt.Lv = Lpval(xtraj_opt.t, xtraj_opt.x(:, 1), xtraj_opt.x(:, 2));
    end

        
    figure(4)        
    clf
    
    %value function 
    subplot(3,1,1)
    
    hold on        
    for i = 1:Nsample
        if i == 1
            plot(xtraj{i}.t(xtraj{i}.x_valid), xtraj{i}.v(xtraj{i}.x_valid)  + obj, 'c','DisplayName', 'Trajectories')
        else
            plot(xtraj{i}.t(xtraj{i}.x_valid), xtraj{i}.v(xtraj{i}.x_valid) + obj, 'c', 'HandleVisibility','off')
        end  
    end

    if rank0 == 1 && rankp == 1
        plot(xtraj_opt.t(xtraj_opt.x_valid), xtraj_opt.v(xtraj_opt.x_valid) + obj, 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
    end
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')

    hold off
    title('Safety Function along Trajectories')
    xlabel('time')
    ylabel('v(x,t) - \gamma')
    legend('location', 'northwest')
    
    subplot(3,1,2)

    
    hold on        
    for i = 1:Nsample
        if i == 1
            plot(xtraj{i}.t(xtraj{i}.x_valid), xtraj{i}.Lv(xtraj{i}.x_valid), 'c','DisplayName', 'Trajectories')
        else
            plot(xtraj{i}.t(xtraj{i}.x_valid), xtraj{i}.Lv(xtraj{i}.x_valid), 'c', 'HandleVisibility','off')
        end  
    end

    if rank0 == 1 && rankp == 1
        plot(xtraj_opt.t(xtraj_opt.x_valid), xtraj_opt.Lv(xtraj_opt.x_valid), 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
    end
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')

    hold off
    title('Change in Safety Function along Trajectories')
    xlabel('time')
    ylabel('L_f v(x,t)')
    legend('location', 'northwest')
    
    subplot(3,1,3)
    
    hold on
    for i = 1:Nsample
        p_comp = -xtraj{i}.cost - xtraj{i}.v;
        
        if i == 1
            plot(xtraj{i}.t(xtraj{i}.x_valid), p_comp(xtraj{i}.x_valid), 'c','DisplayName', 'Trajectories')
        else
            plot(xtraj{i}.t(xtraj{i}.x_valid), p_comp(xtraj{i}.x_valid), 'c', 'HandleVisibility','off')
        end  
    end

    if rank0 == 1 && rankp == 1
        p_comp = -xtraj{i}.cost - xtraj{i}.v - obj;
        plot(xtraj_opt.t(xtraj_opt.x_valid), p_comp(xtraj{i}.x_valid), 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
    end
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')
    hold off
    
    
    title('Comparision with Objective Function')
    xlabel('time')
    ylabel('-cost(x) - v(x,t)')
    legend('location', 'northwest')
    
    
    
end


