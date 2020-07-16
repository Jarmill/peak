%% Parameters for solver    
clear; order = 4; d = 2*order;

mset clear; warning('off','YALMIP:strict')
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));

%% Problem Definition
mpol('x', 2);  %occupation measure    
mpol('xp', 2); %peak measure
DECLARE_MEASURES = 1; %setting this to 0 slows down solving. why?
if DECLARE_MEASURES
    mu = meas(x); mup = meas(xp); 
end

%Moment Constraints
v  = mmon(x, d);    vp = mmon(xp, d); %monomial test functions

powers = genPowGlopti(2, d); %initial measure is a dirac delta, explicit moments
y0 = prod([1.5 0].^powers, 2); yp = mom(vp); %initial set at [1.5, 0]
mom_con = (mom(diff(v, x)*[x(2); -x(1) + (1/3).* x(1).^3 - x(2)]) + (y0 - yp) == 0); %dynamics

%Support Constraints
R = 5;    %radius to contain dynamics
X  = (x'*x <= R^2);     Xp = (xp'*xp <= R^2);
supp_con = [X, Xp];

objective = min(xp(2));

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