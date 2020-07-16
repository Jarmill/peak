%initial set
C0 = [1.5; 0];
R0 = 0.4;

%prajna and rantzer flow
fv = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

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

%moment constraint 
mom_con = [ mom(diff(v, x)*fv(0, x)) + (y0 - yp) == 0; mass(mu0)==1];

%Support Constraints
X  = (x'*x <= R^2);
Xp = (xp'*xp <= R^2);
X0 = ((x0(1)-C0(1))^2 + (x0(2)-C0(2))^2 <= R0^2);

supp_con = [X0, X, Xp];
cost = -xp(2);

objective = max(cost);

%% Solve LMI and recover solution

%Input LMI moment problem
P = msdp(objective, ...
    mom_con, supp_con);
[status,obj,m,dual_rec]= msol(P);  

scale(x, R); scale(x0, R); scale(xp, R);
P_scale = msdp(objective, ...
    mom_con, supp_con);
[status_scale,obj_scale,m_scale,dual_rec_scale]= msol(P_scale);  
