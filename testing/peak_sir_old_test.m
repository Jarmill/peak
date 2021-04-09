% Peak Estimation of SIR system

%Author: Jared Miller 6/25/20
SOLVE = 1;
PLOT = 1;

%C0 = [1; 1];
%C0 = [-0.5; -0.5];
%R0 = 0.5;
%mu_van = 1;

beta = 0.4;
gamma = 0.035 + 0.005;
I_max = 0.1;


%test a sample trajectory
fv = @(t, x) [-beta*x(1)*x(2); beta*x(1)*x(2) - gamma*x(2)];

[tv, xv] = ode45(fv, [0, 20], [1-I_max, I_max]);

cost_max_test = max(xv(:, 2));

if SOLVE
    mset clear
    mset('yalmip',true);
    mset(sdpsettings('solver', 'mosek'));

    d0 = 5;  %degree of relaxation
%     d0 = 2;
    %d0 = input('order of relaxation ='); d = 2*d0;
    %d0 = 10;
    d = 2*d0;
    
    T = 40;   %final time
    R = 4;    %radius to contain dynamics
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
    %f = T*[x(2); mu_van*(1-x(1)^2)*x(2) - x(1)];
    f = T * [-beta*x(1)*x(2); beta*x(1)*x(2) - gamma*x(2)];
    
    %Liouville Equation 
    Ay = mom(diff(v, t) + diff(v, x)*f); 
    Liou = Ay + (y0 - yp);

    mom_con = [Liou == 0; mass(mu0)==1];


    %Support Constraints
    %X  = (x'*x <= R^2);
    %Xp = (xp'*xp <= R^2);
    %X0 = ((x0(1)-C0(1))^2 + (x0(2) - C0(2))^2 <= R0^2);
    
    X = [x(1) >= 0, x(2) >= 0, x(1) + x(2) <= 1];
    Xp = [xp(1) >= 0, xp(2) >= 0, xp(1) + xp(2) <= 1];
    X0 = [x0(1) >= 0, x0(1) + x0(2) <= 1, x0(2) >= 0, x0(2) <= I_max];

    %bounds on trajectory
    supp_con = [X0, X, Xp, ...
        t*(1-t) >= 0, tp*(1 - tp) >= 0, t0 == 0];

    %function to optimize
    cost = xp(2);    
    objective = max(cost);

    %Input LMI moment problem
    P = msdp(objective, ...
        mom_con, supp_con);

    %solve LMIP moment problem
    [status, obj] = msol(P);    
    
    %extract solutions    
    upper_lim = obj;
    M0 = double(mmat(mu0));
    Mp = double(mmat(mup));
    
    %1, t, x1, x2
    M0_1 = M0(1:4, 1:4);
    Mp_1 = Mp(1:4, 1:4);
    
    rankp = rank(Mp_1, 1e-3);
    rank0 = rank(M0_1, 1e-3);
    
    t0_out = double(mom(t0));
    x0_out = double(mom(x0));
    
    tp_out = T*double(mom(tp));
    xp_out = double(mom(xp));    
    
end

if PLOT        
    %sample from X0 
    Nsample = 80;
    
    
    rng(33, 'twister')
    X0_pts = [0, 1, 1-I_max, 0, 0;
          0, 0, I_max, I_max, 0];
    
    
    Xsample = trap_sample(Nsample, I_max);    

                   
    
        
    figure(1)
    clf
    subplot(1,3,1)
    hold on
    plot(X0_pts(1, :), X0_pts(2, :), 'k', 'Linewidth', 3)
    plot([0, 1], [upper_lim, upper_lim], 'r--', 'Linewidth', 3)
    
    xtraj = cell(Nsample, 1);
    for i = 1:Nsample
        xtraj{i} = struct;
        [xtraj{i}.t, xtraj{i}.x] = ode45(fv, [0, T], [Xsample(:, i)]);  
        
        if i == 1
            plot(xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c')
        else
            plot(xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c', 'HandleVisibility','off')
        end
    end 
    
    scatter(Xsample(1, :), Xsample(2, :), 'k')
    
    if rank0 == 1 && rankp == 1
        %global optimum was found
        
        xtraj_opt = struct;
        [xtraj_opt.t, xtraj_opt.x] = ode45(fv, [0, T], x0_out);  
        
        plot(xtraj_opt.x(:, 1), xtraj_opt.x(:, 2), 'b', 'LineWidth', 2)
        
        MS = 200;
        
        scatter(xp_out(1), xp_out(2), MS, '*b', 'HandleVisibility','off', 'LineWidth', 2)
        scatter(x0_out(1), x0_out(2), MS, 'ob', 'HandleVisibility','off', 'LineWidth', 2)
        
        legend({'Initial Set', 'Upper Limit', 'Trajectories', 'Peak Traj.'}, 'location', 'northeast')
        title(['Peak Infection Rate = ', num2str(obj, 3), ' at order = ', num2str(d0)])
    else
        %global optimum not certified
        legend({'Initial Set', 'Upper Limit', 'Trajectories'}, 'location', 'northeast')
        title(['Upper Limit for Infection Rate = ', num2str(obj, 3), ' order = ', num2str(d0)])
    end
    
    
    
    
    
    xlim([0, 1])
    ylim([0, 1])
    hold off
    axis square
    xlabel('Susceptible')
    ylabel('Infected')
    
    subplot(1,3,[2, 3])    
    hold on
    
    plot([0, 0], [0, I_max], 'k', 'LineWidth', 3);
    plot([0, T], [upper_lim, upper_lim], 'r--', 'Linewidth', 3)
    
    for i = 1:Nsample        
        if i == 1
            plot(xtraj{i}.t, xtraj{i}.x(:, 2), 'c')
        else
            plot(xtraj{i}.t, xtraj{i}.x(:, 2), 'c', 'HandleVisibility','off')
        end
    end 
    
    if rank0 == 1 && rankp == 1
        %global optimum was found
        %hlines_p = streamline(x, y, xdot, ydot, x0_out(1), x0_out(2));
        %set(hlines_p, 'LineWidth', 2)
        
        
        plot(xtraj_opt.t, xtraj_opt.x(:, 2), 'b', 'LineWidth', 2)
        
        
        scatter(tp_out, xp_out(2), MS, '*b', 'HandleVisibility','off', 'LineWidth', 2)
        scatter(0, x0_out(2), MS, 'ob', 'HandleVisibility','off', 'LineWidth', 2)
        
        legend({'Initial Set', 'Upper Limit', 'Trajectories', 'Peak Traj.'}, 'location', 'northeast')
        
        title(['Peak Infection Rate = ', num2str(obj, 3), ' at time = ', num2str(tp_out, 3)])
    else
        %global optimum not certified
        legend({'Initial Set', 'Upper Limit', 'Trajectories'}, 'location', 'northeast')        
        title(['Upper Limit for Infection Rate = ', num2str(obj, 3)])
    end
    
    %axis square
    xlabel('Time')
    ylabel('Infection Rate ')
    pbaspect([2.25 1 1])

    
end

function [X] = circle_sample(N)
%CIRCLE_SAMPLE Uniformly Sample N points from the unit circle

theta = 2*pi*rand(N, 1);
trig = [cos(theta) sin(theta)];
r = sqrt(rand(N, 1));

X = r.*trig;
end

function [X] = trap_sample(N, I_max)
% sample from trapezoid [0<=x<=1], [0<=y<=I_max]

x = rand(15*N, 1);
y = rand(15*N, 1) * (I_max);

rej = (x+y) > 1;

x(rej) = [];
y(rej) = [];

X = [x(1:N), y(1:N)]';


end
