function [peak_val, out] = peak_impulse_siso(A, B, C, order, rank_tol, Tmax)
%PEAK_IMPULSE_SISO  Find the maximum value of |Cx| for the impulse response
%of a SISO linear system that switches between dynamics x'=Ai x+Bu, y=Cx.
%arbirarily for a set of matrices Ai (simple case where there is only one
%Ai, standard linear dynamics
%Input:
%   A:          Cell of plausible linear dynamics (or just one matrix)
%   B:          Impulse response -> Initial condition
%   C:          Output to optimize, max |Cx|
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

if ~iscell(A)
    A = {A};
end

%Parameters of relaxation
n = size(A{1}, 1);
nsys = length(A);
d = order*2;

%I know this is invali
%maxsvd = max(cellfun(@(Ai) svds(Ai, 1), A));

%Gloptipoly Options
mset clear; warning('off','YALMIP:strict')
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));


%% Set up the measures
mpol('x', n, nsys);  
%mpol('x', n, 1);  
mpol('xp', n);  


%assume Tmax == Inf, makes life easy

%if Tmax == Inf
%    var = x;
%    varp = xp;
%else
%    mpol('t', 1, nsys);
%    mpol('tp', 1);
%    var = [t; x];
%    varp = [tp; xp];
%end

%measures


%occupation measure for each subsystem
mu = cell(nsys, 1);
v = cell(nsys, 1);

%Liouville
%Ayi = cell(nsys, 1);
Ay = 0;

%scaling and support
R = 100; %something large, since joint spectral radius plays havoc
         %as R reduces, solution grows more accurate
X = [];
for i = 1:nsys   
    xcurr = x(:, i);
    %xcurr = x;
    mu{i} = meas(xcurr);
    v{i}  = mmon(xcurr, d);
    
    %Liouville, mu = sum_mi, L'mu = sum_i Lf' mu_i
    %Ayi{i} = mom(diff(v, x)*A*x);
    Ay = Ay + mom(diff(v{i}, xcurr)*A{i}*xcurr);
    
    %support
    X = [X, xcurr'*xcurr <= R];
end
%mu  = meas(var);   


%peak measure
%varp = xp;
%mup = meas(varp); 
mup = meas(xp);
vp = mmon(xp, d);
yp = mom(vp);
%support
Xp = (xp'*xp <= R^2);

% 
% if Tmax ~= Inf
%     Ay = mom(diff(v, t)) + Ay*Tmax;
% end

%initial measure
% if Tmax == Inf
    powers = genPowGlopti(n, d);
    y0 = prod((B').^powers, 2);
% else
%     %starts at time 0
%    powers = genPowGlopti(n+1, d);
%     y0 = prod(([0; B]').^powers, 2);
% end

%% Gloptipoly constraints
Liou = Ay + (y0 - yp);

mom_con = (Liou == 0);

%support
%maximum norm, approximate with magnification (max svd) of A.
%R = 10*maxsvd*norm(B)^2;


%X  = (x'*x <= R^2);
% Xp = (xp'*xp <= R^2);
% 
% if Tmax ~= Inf
%     X  = [X, t*(1-t) >= 0];
%     Xp = [Xp, tp*(1-tp) >= 0];
% end
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

% if Tmax == Inf
    vv = conj(monolistYalToGlop(xc, d));
    p = dual_rec'*vv;
    Lp = cell(nsys, 1);
    for i = 1:nsys
        Lp{i} = jacobian(p, xc)*A{i}*xc;
    end
    
    
% else
%     syms tc;
%     
%     vv = conj(monolistYalToGlop([T*tc, xc], d));
%     p = dual_rec'*vv;
%     Lp = diff(p, tc) + jacobian(p, xc)*A*xc;
% end


pval = matlabFunction(p);
if nsys == 1
    Lpval = matlabFunction(Lp{1});
else
    Lpval = cellfun(@(Lpf) matlabFunction(Lpf), Lp, 'UniformOutput', false);
end

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

