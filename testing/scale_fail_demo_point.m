%% Parameters for solver    
mset clear; warning('off','YALMIP:strict')
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));

order = 5; d = 2*order;

%% Problem Definition
%measures
mpol('x', 2);   %occupation measure    
mpol('xp', 2); %peak measure

%Support Constraints
R = 5;    %radius to contain dynamics
X  = (x'*x <= R^2);
Xp = (xp'*xp <= R^2);

supp_con = [X, Xp];

cost = -xp(2); %minimum height of trajectory
objective = max(cost);

%Measures, test functions (monomials)
v  = mmon(x, d);    vp = mmon(xp, d);

%initial set at [1.5, 0]
powers = genPowGlopti(2, d);
y0 = prod([1.5 0].^powers, 2);  yp = mom(vp);

%moment constraint, Liouville and dynamics
mom_con = (mom(diff(v, x)*[x(2); -x(1) + (1/3).* x(1).^3 - x(2)]) + (y0 - yp) == 0);

%% Solve LMI and recover solution

%Input LMI moment problem
P = msdp(objective, ...
    mom_con, supp_con);
[status,obj,m,dual_rec]= msol(P);  


%with scaling
scale(x, R);  scale(xp, R);
P_scale = msdp(objective, ...
    mom_con, supp_con);
[status_scale,obj_scale,m_scale,dual_rec_scale]= msol(P_scale);  

%scale difference
scale_diff = obj_scale - obj;