% Hybrid Peak Estimation: height of bouncing ball
% based on 'ballode' demo
% Author: Jared Miller 1/23/21
SOLVE = 1;
PLOT = 1;

%use prior knowledge that vertical velocity=0 at peak?
zero_v_at_peak = 1;

T = 10;
g = 1;

%reset map
kappa = 0.9; %speed preserved after bounce
Nr = 10; %number of bounces permitted (mass of measure)

%throw ball against the ground, observe maximum height
x0_orig = [1; -2];

if SOLVE
    mset clear
    mset('yalmip',true);
    mset(sdpsettings('solver', 'mosek'));

    %d = 2*2;  %degree of relaxation
    d0 = 3;
    %d0 = input('order of relaxation ='); d = 2*d0;
    %d0 = 10;
    d = 2*d0;
    
    
    %variables and measures
    mpol('t', 1);  mpol('x', 2);  mu  = meas([t; x]);   %occupation measure
    mpol('t0', 1); mpol('x0', 2); mu0 = meas([t0; x0]); %initial measure
    mpol('tp', 1); mpol('xp', 2); mup = meas([tp; xp]); %peak measure
    mpol('tr', 1); mpol('xr', 2); mur = meas([tr; xr]); %reset (counting) measure

    %test functions (monomials)
    v  = mmon([t; x], d);
    v0 = mmon([t0; x0], d);
    vp = mmon([tp; xp], d);
    vr = mmon([tr; xr], d);

    %unknown moments of initial measure
    y0 = mom(v0);
    yp = mom(vp);

    f = T * [x(2); -g];
   
    %kappa: jacobian factor. May be 1/kappa? test out_traj.
%     R = [xr(1); -kappa*xr(2)];
%     DR = diff(R, xr);
%     detR = det(DR); %absolute value of jacobian?
%     detR = kappa;
    %no jacobian needed
    reset_push = subs(vr, xr, [xr(1); -kappa*xr(2)]);
    %Liouville Equation 
    Ay = mom(diff(v, t) + diff(v, x)*f); 
    Ar = mom(reset_push - vr);
    
    Liou = Ar + Ay + (y0 - yp);

    mom_con = [Liou == 0; mass(mu0)==1; mass(mur)<=Nr];


    %Support Constraints
    X  = [x(1) >= 0, x(2)^2 <= 9];
    if zero_v_at_peak
        %prior information
        %at peak of ball motion, there is no vertical velocity
        Xp = [xp(1) >= 0, xp(2) == 0]; %may be able to do xp(2) == 0   
    else
        Xp = [xp(1) >= 0, xp(2)^2 <= 9]; %may be able to do xp(2) == 0   
    end
    X0 = [x0 == x0_orig]';

    Xr = [xr(1) == 0, xr(2)^2 <= 9];

    %bounds on trajectory
    supp_con = [X0, X, Xp, Xr, ...
        t*(1-t) >= 0, tp*(1 - tp) >= 0, tr*(1 - tr) >= 0, t0 == 0];

    %function to optimize
    cost = xp(1);    %height of ball
    objective = max(cost);

    %Input LMI moment problem
    P = msdp(objective, ...
        mom_con, supp_con);

    %solve LMIP moment problem
    [status, obj] = msol(P);    
    
    %extract solutions    
    upper_lim = obj
    M0 = double(mmat(mu0));
    Mp = double(mmat(mup));
    Mocc = double(mmat(mu));
    Mr = double(mmat(mur));
    %1, t, x1, x2
    M0_1 = M0(1:4, 1:4);
    Mp_1 = Mp(1:4, 1:4);
    
    rankp = rank(Mp_1, 1e-3);
    rank0 = rank(M0_1, 1e-3);
    
    t0_out_traj = double(mom(t0));
    x0_out_traj = double(mom(x0));
    
    tp_out_traj = T*double(mom(tp));
    xp_out_traj = double(mom(xp));    
    
end

if PLOT
    out_traj = bouncing_ball(x0_orig, T, kappa);
    
    fig = figure(1);
    clf
    subplot(3, 1, [1, 2])
%     box on
    hold on;
    plot(T*out_traj.t,out_traj.y(:,1),'-')
    scatter(T*out_traj.te,out_traj.ye(:,1),150,'k','o')
    plot([0 T], [obj obj], 'r--', 'LineWidth', 2)
    if rankp == 1
        scatter(tp_out_traj, xp_out_traj(1), 150, 'r', '*')
    end
    xlim([0, T]);
    ylim([0, 3.5]);
    xlabel('time (s)');
    ylabel('height (m)');
    title(['Ball height bound = ', num2str(obj, 3), ', degree = ', num2str(d)], 'Fontsize', 14);
    hold off
    
    subplot(3, 1, 3)
    hold on
    plot(T*out_traj.t, out_traj.y(:, 2))
    scatter(T*out_traj.te,out_traj.ye(:,2),150,'k','o')
    ylabel('velocity (m/s)')
    xlabel('time (s)')
    title('Ball velocity', 'Fontsize', 14)
    if rankp == 1
        scatter(tp_out_traj, xp_out_traj(2), 150, 'r', '*')
    end
    hold off
end