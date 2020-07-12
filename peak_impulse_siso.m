function [peak_val, out] = peak_impulse_siso(A, B, C, order, rank_tol, Tmax)
%PEAK_IMPULSE_SISO  Find the maximum value of |Cx| for the impulse response
%of a SISO linear system with dynamics x'=Ax+Bu, y=Cx.
%Input:
%   A,B,C:      Linear System
%   order:      Order of moment relaxation (degree=2*order)
%   ranktol:    rank tolerance of moment matrixof 
%   Tmax:       Maximum time. By default is Infinite (time-independent).
%               Tmax ~= Inf will recover optimal time at optimality, and
%               allow for time-varying safety contours
%Ouputs:
%   peak_val:   Maximum absolute value of impulse response
%   out:        Data structure that contains the peak response information   

if nargin < 6
    Tmax = Inf;    
end

if nargin < 5
    rank_tol = 1e-3;
end

if nargin < 4
    order = 2;
end

%Parameters of relaxation
n = size(A, 1);
d = order*2;

%Gloptipoly Options
mset clear; warning('off','YALMIP:strict')
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));


%% Set up the measures
mpol('x', n);  
mpol('xp', n);  


if Tmax == Inf
    var = x;
    varp = xp;
else
    mpol('t', 1);
    mpol('tp', 1);
    var = [t; x];
    varp = [tp; xp];
end

%measures
mu  = meas(var);   %occupation measure    
mup = meas(varp); %peak measure

v  = mmon(x, d);
vp = mmon(xp, d);
yp = mom(vp);

%Liouville

Ay = mom(diff(v, x)*A*x);

if Tmax ~= Inf
    Ay = mom(diff(v, t)) + Ay*Tmax;
end

%initial measure
if Tmax == Inf
    powers = genPowGlopti(n, d);
    y0 = prod((B').^powers, 2);
else
    %starts at time 0
    powers = genPowGlopti(n+1, d);
    y0 = prod(([0; B]').^powers, 2);
end

%% Gloptipoly constraints
Liou = Ay + (y0 - yp);

mom_con = (Liou == 0);

%support
%maximum norm, approximate with magnification (max svd) of A.
R = 10*svds(A, 1)*norm(B)^2;

X  = (x'*x <= R^2);
Xp = (xp'*xp <= R^2);

if Tmax ~= Inf
    X  = [X, t*(1-t) >= 0];
    Xp = [Xp, tp*(1-tp) >= 0];
end
supp_con = [X, Xp];

cost = (C*xp)^2;

objective = max(cost);


%% Solve the Problem
cost = (C*xp)^2;

objective = max(cost);

P = msdp(objective, ...
    mom_con, supp_con);



%solve LMI moment problem    
[status,obj,m,dual_rec]= msol(P);

%Extract Moments (Primal)
peak_val = sqrt(obj);

Mp = double(mmat(mup));
xp_out = double(mom(xp));

if Tmax == Inf
    Mp_1 = Mp(1:3, 1:3);
    tp_out = Inf;
else
    Mp_1 = Mp(1:4, 1:4);
    tp_out = double(mom(tp));
    
end
    
rankp = rank(Mp_1, rank_tol);

%Extract Polynomials (Dual)

syms xc [n, 1];

if Tmax == Inf
    vv = conj(monolistYalToGlop(xc, d));
    p = dual_rec'*vv;
    Lp = jacobian(p, xc)*A*xc;
else
    syms tc;
    
    vv = conj(monolistYalToGlop([T*tc, xc], d));
    p = dual_rec'*vv;
    Lp = diff(p, tc) + jacobian(p, xc)*A*xc;
end


pval = matlabFunction(p);
Lpval = matlabFunction(Lp);

%% package together outputs
out = struct;
out.peak_val = peak_val;
out.obj = obj;
out.xp = xp_out;
out.tp = tp_out;
out.Mp_1 = Mp_1;
out.Mp = Mp;
out.rankp = rankp;

%level set: out.p + out.obj = 0

%polynomials
out.p = p;
out.pval = pval;
out.Lp = Lp;
out.Lpval = Lpval;
out.cost = @(x) (x*C').^2;

end

