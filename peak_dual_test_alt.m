%solve the flow test program in dual form
%hopefully this works

SOLVE = 0;
PLOT = 1;
ISOSURFACE = 1;

%initial set
C0 = [1.5; 0];
R0 = 0.5;

%unsafe set
Cu = [-1; -1];
Ru = 0.4;


%minimize the vertical coordinate (1)
%or safety margin of circle (0)
LINE_COST = 0;

%mesh density
%MD = 200;
MD = 80;

%final time
T = 6;
%dynamics
fv = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

if SOLVE
    %options and degrees
    opts = sdpsettings('solver', 'mosek');
    order =4;
    d = order * 2;

    %variables in polynomials
    x = sdpvar(2, 1);
    t = sdpvar(1, 1);

    %define polynomial v and parameter gamma
    [v, cv] = polynomial([t; x], d);
    gamma = sdpvar(1, 1);

    %constraint set
    R = 5;
    
    %dynamics
    %f = T*[x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
    f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

    Lv = jacobian(v, t) + jacobian(v, x) * f;

    % unsafe region level
    %p =  (x(1)+1)^2 + (x(2)+1)^2 - 0.16;
    if LINE_COST
        p = -x(2);   
    else
        p =  Ru^2 -(x(1)-Cu(1))^2 - (x(2)-Cu(2))^2;
    end
    %Initial set
    gX0 = R0^2 - (x(1)-C0(1))^2 - (x(2) - C0(2))^2;
    gX = R^2 - x(1)^2 - x(2)^2;
    %gt = t*(1-t);
    gt = t*(T-t);

    %Putinar Multipliers
    [sx0, cx0] = polynomial(x, d-2);
    [stf, ctf] = polynomial([t;x], d-2);
    [sxf, cxf] = polynomial([t;x], d-2);
    [stp, ctp] = polynomial([t;x], d-2);
    [sxp, cxp] = polynomial([t;x], d-2);

    obj = -gamma;

    %dual program:

    %cost:   max gamma
    %mu0:    v(0, x) >= gamma      x in X0
    %mu:     Lv(t,x) >= 0          t in [0,T], x in X
    %mup:    p(t,x)  >= v(t,x)     t in [0,T], x in X

    %initial conditions
    cons = [sos(replace(v, t, 0) - gamma - gX0*sx0), sos(sx0)];
    %dynamics
    cons = [cons, sos(Lv - gt*stf - gX*sxf), sos(stf), sos(sxf)];
    %peak
    %cons = [cons, sos(p - v - gt*stp - gX*sxp), sos(stp), sos(sxp)];
    cons = [cons, sos(-p - v - gt*stp - gX*sxp), sos(stp), sos(sxp)];

    solvesos(cons, obj, opts, [cv; cx0; ctf; cxf; ctp; cxp; gamma]);

    obj_val = value(obj);
    gamma_val = value(gamma);
    radius_out = sqrt(-value(obj) + Ru^2);
end

if PLOT
    m = 3.25;
    Nsample = 60;
    %initial set        
    Ntheta = 100;
    theta = linspace(0, 2*pi, Ntheta);
    circ = [cos(theta); sin(theta)];
    
    X0 = C0 + circ*R0;
        
    
    %sample from X0 
    rng(33, 'twister')

    Xsample = C0 + circle_sample(Nsample)' * R0;
    
    Xu = Cu + circ*Ru;            
    Xmargin = Cu + circ*radius_out;
    
    
    figure(3)
    clf
    %plot the phase plane
    subplot(1, 3, 1)    
    hold on
    plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3)
    
    if LINE_COST
        plot([-m, m], gamma_val*[1, 1], 'r--', 'Linewidth', 3)
    else
        patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
        plot(Xmargin(1, :), Xmargin(2, :), 'r--', 'Linewidth', 3)
    end
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
    
    xlim([-m, m])
    ylim([-m, m])
    axis square
    xlabel('x_1')
    ylabel('x_2')
    
    if LINE_COST
        legend({'Initial Set', 'Escape Margin', 'Trajectories'}, 'location', 'northwest')
        title(['Escape Margin for Trajectories = ', num2str(gamma_val, 2), ' order = ', num2str(order)])
    else
    legend({'Initial Set', 'Unsafe Set', 'Escape Margin', 'Trajectories'}, 'location', 'northwest')
        title(['Escape Margin for Trajectories = ', num2str(gamma_val, 2)])
    end
    
    %plot trajectories in time
    subplot(1,3,[2,3])    
    hold on
    
    %initial set
    plot3(zeros(Ntheta, 1), X0(1, :), X0(2, :), 'k', 'Linewidth', 3)
    
    %lower bound to maximum (escape margin)
    if LINE_COST                
        patch('XData',[0, 0, T, T], 'YData',[-m, m, m, -m], 'ZData', gamma_val*[1,1,1,1], 'FaceColor', 'r', 'FaceAlpha', 0.3)
    else
        [xcyl, ycyl, zcyl] = cylinder(1, 50);
        surf(T*zcyl, Cu(1) + xcyl*Ru, Cu(2) + ycyl*Ru, 'FaceColor', 'r', 'EdgeColor', 'None')
        surf(T*zcyl, Cu(1) + xcyl*radius_out, Cu(2) + ycyl*radius_out, ...
            'FaceColor', 'r', 'EdgeColor', 'None', 'Facealpha', 0.3)
    end
    
    %contour in SOS program (escape contour)
    if ISOSURFACE
        syms tc xc yc;
        vv = monolist([tc; xc; yc], d);
        p = value(cv)'*vv - value(gamma);

        
        fi = fimplicit3(p, [0, T, -m, m, -m, m], 'MeshDensity',MD, ...
            'EdgeColor', 'None','FaceColor', 'k', 'FaceAlpha', 0.3, ...
            'DisplayName', 'Escape Contour');

    end
    
    for i = 1:Nsample        
        if i == 1
            plot3(xtraj{i}.t, xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c')
        else
            plot3(xtraj{i}.t, xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c', 'HandleVisibility','off')
        end        
    end  
    
    if ISOSURFACE
        legend({'Initial Set', 'Escape Margin', 'Escape Contour', 'Trajectories'}, 'location', 'northwest')
    else
        legend({'Initial Set', 'Escape Margin', 'Trajectories'}, 'location', 'northwest')
    end
    
    title(['Escape Margin for Trajectories = ', num2str(gamma_val, 2), ' order = ', num2str(order)])
    
    xlim([0, T])
    ylim([-m, m])
    zlim([-m, m])
    xlabel('time')
    ylabel('x_1')
    zlabel('x_2')
    pbaspect([2.25 1 1])
    
    
    
end


function [X] = circle_sample(N)
%CIRCLE_SAMPLE Uniformly Sample N points from the unit circle

theta = 2*pi*rand(N, 1);
trig = [cos(theta) sin(theta)];
r = sqrt(rand(N, 1));

X = r.*trig;
end

%plots?
%work on that next