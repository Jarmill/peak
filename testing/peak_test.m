% Try the peak estimation on a vector field
% Dynamics from 'On Analysis and Synthesis of Safe Control Laws'
% by Anders Rantzer and Stephen Prajna

%Author: Jared Miller 6/22/20
SOLVE = 1;
PLOT = 1;

T = 10;   %final time

%initial set
C0 = [1.5; 0];
%R0 = 0.5;
R0 = 1;

%unsafe set
Cu = [-1; -1];
Ru = 0.4;

%dynamics
fv = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

if SOLVE
    mset clear; warning('off','YALMIP:strict')
    mset('yalmip',true);
    mset(sdpsettings('solver', 'mosek'));

    %d = 2*2;  %degree of relaxation
    %d0 = input('order of relaxation ='); d = 2*d0;
    d0 = 2;
    d = 2*d0;
        
    R = 5;    %radius to contain dynamics
    %dynamics are time-independent

    %variables and measures
    mpol('t', 1);  mpol('x', 2);  mu  = meas([t; x]); %occupation measure
    mpol('t0', 1); mpol('x0', 2); mu0 = meas([t0; x0]);
    mpol('tp', 1); mpol('xp', 2); mup = meas([tp; xp]);

    %test functions (monomials)
    v  = mmon([t; x], d);
    v0 = mmon([t0; x0], d);
    vp = mmon([tp; xp], d);

    %unknown moments of initial measure
    y0 = mom(v0);
    yp = mom(vp);

    %dynamics and cost
    %scale all dynamics to [0, 1] to improve time scaling
    f = T*[x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
    %Liouville Equation 
    Ay = mom(diff(v, t) + diff(v, x)*f); 
    Liou = Ay + (y0 - yp);

    mom_con = [Liou == 0; mass(mu0)==1];


    %Support Constraints
    X  = [x'*x <= R^2];
    Xp = [xp'*xp <= R^2];
    X0 = [(x0(1)-C0(1))^2 + (x0(2)-C0(2))^2 <= R0^2];

    %bounds on trajectory
    %supp_con = [X0, X, Xp, ...
        %t*(1-t) >= 0, tp*(1 - tp) >= 0, t0 == 0];
    
    supp_con = [X0, t0 == 0, X, t*(1-t) >= 0, Xp, tp*(1 - tp) >= 0];
        
    %function to optimize
    cost = Ru^2 - (xp(1)-Cu(1))^2 - (xp(2)-Cu(2))^2;
    %cost = -xp(2);
    objective = max(cost);

    %Input LMI moment problem
    P = msdp(objective, ...
        mom_con, supp_con);

    %solve LMIP moment problem
    %[status, obj] = msol(P);
    [status,obj,m,dual]= msol(P);

    
    %extract solutions
    radius_out = sqrt(-obj + Ru^2);
    M0 = double(mmat(mu0));
    Mp = double(mmat(mup));
    
    %1, t, x1, x2
    M0_1 = M0(1:4, 1:4);
    Mp_1 = Mp(1:4, 1:4);
    
    rank0 = rank(Mp_1, 1e-4);
    rankp = rank(M0_1, 1e-4);
    
    t0_out = double(mom(t0));
    x0_out = double(mom(x0));
    
    tp_out = T*double(mom(tp));
    xp_out = double(mom(xp));
    %radius_out = sqrt(obj);
    
    [model.A, model.b, model.c, model.K, model.b0, model.s] = msedumi(P);
    
    
    
end

if PLOT
    m = 3;
    %N = 20;
    Nsample = 40;
    %[x, y] = meshgrid(linspace(-m, m, N));

    %xdot = y;
    %ydot = -x + (1/3).* x.^3 - y;

    %quiver(x, y, xdot, ydot)

    %initial and unsafe sets
    theta = linspace(0, 2*pi, 100);
    circ = [cos(theta); sin(theta)];
    
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
    hold on
    
    plot(X0_pts(1, :), X0_pts(2, :), 'k', 'Linewidth', 3)
    patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
    plot(Xmargin(1, :), Xmargin(2, :), 'r--', 'Linewidth', 3)
    %hlines_c0 = streamline(x, y, xdot, ydot, C0(1), C0(2));
    %hlines = streamline(x, y, xdot, ydot, Xsample(1, :), Xsample(2, :));
    
    %set(hlines_c0, 'Color', 'c')
    %set(hlines, 'Color', 'c', 'HandleVisibility','off')
    
    xtraj = cell(Nsample, 1);
    for i = 1:Nsample
        xtraj{i} = struct;
        [xtraj{i}.t, xtraj{i}.x] = ode45(fv, [0, T], Xsample(:, i));  
        
        if i == 1
            plot(xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c')
        else
            plot(xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c', 'HandleVisibility','off')
        end
    end 
    
    
    if rank0 == 1 && rankp == 1
        %global optimum was found
        
        xtraj_opt = struct;
        [xtraj_opt.t, xtraj_opt.x] = ode45(fv, [0, T], x0_out);  
        
        plot(xtraj_opt.x(:, 1), xtraj_opt.x(:, 2), 'b', 'LineWidth', 2)
        
        
        MS = 200;
        
        scatter(xp_out(1), xp_out(2), MS, '*b', 'HandleVisibility','off', 'LineWidth', 2)
        scatter(x0_out(1), x0_out(2), MS, 'ob', 'HandleVisibility','off', 'LineWidth', 2)
        
        legend({'Initial Set', 'Unsafe Set', 'Safety Margin', 'Trajectories', 'Peak Traj.'}, 'location', 'northwest')
        title(['Safety Margin for Trajectories = ', num2str(obj, 3), ' at time t=', num2str(tp_out, 3)])
    else
        legend({'Initial Set', 'Unsafe Set', 'Safety Margin', 'Trajectories'}, 'location', 'northwest')
        title(['Safety Margin for Trajectories = ', num2str(obj, 2)])
    end                    

    xlim([-m, m])
    ylim([-m, m])
    
    hold off
    axis square
    xlabel('x')
    ylabel('y')
    
end

function [X] = circle_sample(N)
%CIRCLE_SAMPLE Uniformly Sample N points from the unit circle

theta = 2*pi*rand(N, 1);
trig = [cos(theta) sin(theta)];
r = sqrt(rand(N, 1));

X = r.*trig;
end


