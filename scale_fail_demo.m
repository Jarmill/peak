SOLVE = 1;
PLOT = 1;

SCALE = 1;

%initial set
C0 = [1.5; 0];
R0 = 0.4;

%unsafe set
Cu = [-1; -1];
Ru = 0.4;

%plotting
m_low = -3;
m_high = 3;

%prajna and rantzer flow
fv = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

LINE_COST = 0;

%subspace angle
theta = 3*pi/2; %(equivalent to maximizing -x(2))
%theta = pi; %max -x(1)
%theta = 11*pi/8;
%theta = 3*pi/4;
LINE_COST = 1;


if SOLVE
    %% Parameters for solver    
    mset clear; warning('off','YALMIP:strict')
    mset('yalmip',true);
    mset(sdpsettings('solver', 'mosek'));

    %order = input('order of relaxation ='); d = 2*order;
    order = 5; d = 2*order;
    
    R = 5;    %radius to contain dynamics
    
    
    %measures
    mpol('x0', 2); mu0 = meas(x0); %initial measure
    mpol('x', 2);  mu  = meas(x);   %occupation measure    
    mpol('xp', 2); mup = meas(xp); %peak measure

    %test functions (monomials)

    v0 = mmon(x0, d);
    v  = mmon(x, d);
    vp = mmon(xp, d);

    %unknown moments of initial measure
    y0 = mom(v0);
    yp = mom(vp);

    t = 0;

    f = fv(t, x);

    %Liouville Equation 
    Ay = mom(diff(v, x)*f); 


    if SCALE
        scale(x, R); scale(x0, R); scale(xp, R);
    end

    
        %Liouville Equation
    Liou = Ay + (y0 - yp);

    mom_con = [Liou == 0; mass(mu0)==1];


    %Support Constraints
    X  = (x'*x <= R^2);
    Xp = (xp'*xp <= R^2);
    X0 = ((x0(1)-C0(1))^2 + (x0(2)-C0(2))^2 <= R0^2);
   
        
    supp_con = [X0, X, Xp];
    
    if LINE_COST                
        costv = @(x) x(:, 2)*sin(theta) + x(:, 1)*cos(theta);
        %x * [cos(theta); sin(theta)];
    else
        costv = @(x) Ru^2 - (x(:, 1)-Cu(1)).^2 - (x(:, 2)-Cu(2)).^2;
    end
    
    cost = costv(xp');
    
    objective = max(cost);

    %% Solve LMI and recover solution
    
    %Input LMI moment problem
    P = msdp(objective, ...
        mom_con, supp_con);

    %solve LMI moment problem    
    [status,obj,m,dual_rec]= msol(P);
   
    %extract solutions
    gamma_val = -obj;
    if LINE_COST
        radius_out = sqrt(gamma_val + Ru^2);
    end
    M0 = double(mmat(mu0));
    Mp = double(mmat(mup));
    
    M0_1 = M0(1:3, 1:3);
    Mp_1 = Mp(1:3, 1:3);
    
    rank0 = rank(M0_1, 1e-4);
    rankp = rank(Mp_1, 1e-4);
    
    %space moments
    x0_out = double(mom(x0));        
    xp_out = double(mom(xp));
end