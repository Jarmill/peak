function [out] = peak_estimate_dual(options, order)
%PEAK_ESTIMATE_DUAL  Find the maximum value of a function p(x) along all
%trajectories starting from an initial set X0 along polynomial dynamics
%f(t, x). This is for finite time [0, Tmax], but if f(t, x) is autonomous
%then the problem can be solved for infinite time.
%   Uses the Dual SOS formulation (rename better later)
%   YALMIP is unable to recover moments of optimum solution.
%   Accomodates switching as well.
%
%Input: options structure (peak_options) with fields:
%   var:        Structure of symbolic variables (sdpvar)
%       t:      time (default empty)
%       x:      state
%       w:      parametric uncertainty (default empty)
%
%   Tmax:       Maximum time (if var.t is not empty)
%
%   dynamics:   Structure with fields (f, X)
%       f:      dynamics x' = f(t,x) over the space X.
%                   Each entry is a polynomial                
%       X:      over what space do the dynamics evolve
%                   this allows for switching, in case f and X are cells.
%                   Each entry is an array of polynomials
%
%   obj:        Functions to maximize along trajectories
%       Scalar:     Single function
%       List:       Maximize the minimum of all objectives
%                      Each entry is a polynomial       
%   
%   state_supp: Support set of X
%   state_init: Support set of X0 (initial)
%
%   param:      Support of parametric uncertainty w (if var.w nonempty)
%       
%   rank_tol:   Rank tolerance for moment matrix to be rank-1
%
%Input (separate argument for convenience)
%   order:      Order of relaxation
%
%Output: out structure with fields
%   peak_val:   Value of estimated peak optimum
%   optimal:    Is peak_val a global-optimum with rank-tolerance rank_tol
%
%   x0:         Optimal initial condition
%   x:         Optimal point (obj(x) = peak_val)
%   tp:         Time of optimal point (if time-dependent, otherwise Inf)
%   Mfix:       Moment Matrices of Fixed Distributions
%   Mp:         Moment Matrix of Peak Distribution
%
%   v:          Polynomial from dual program, v(t, x) - gamma >=0 invariant
%   Lv:         Lie derivatives along trajectories of system for v

if nargin < 2
    order = 3;
end
d = 2*order;

%% Process option data structure

%state variable
nx = length(options.var.x);
nvar = nx;

%parameter variables
nw = length(options.var.w);
nvar = nvar + nw;
if nw > 0
    w = options.var.w;
    W = options.param;
else
    w = [];
    W = [];
end
W = fill_constraint(W);


%compact support (artificial)

%scaling of the box
x = options.var.x;

XR = options.R^2 - x'*x;

%deal with boxes and scaling later


%number of switching subsystems
f = options.dynamics.f;
X = options.dynamics.X;

Xall = fill_constraint(options.state_supp);
Xall.ineq = [Xall.ineq; W.ineq; XR];
Xall.eq = [Xall.eq; W.eq];

X0 = fill_constraint(options.state_init);

%X0.ineq = [X0.ineq; W.ineq; XR];
X0.ineq = [X0.ineq; W.ineq];
X0.eq = [X0.eq; W.eq];
%leave risky, make sure that X0 is compact

if iscell(options.dynamics.f)
    nsys = length(options.dynamics.f);
else
    nsys = 1;
    options.dynamics.f = {options.dynamics.f};
    options.dynamics.X = {options.dynamics.X};    
    f = {f};
    X = {X};
end


%fill in equalities and inequalities (ineq >= 0) if not given
for i = 1:nsys    
    X{i} = fill_constraint(X{i});
    
    %add redundant constraint
    X{i}.ineq = [Xall.ineq; X{i}.ineq; W.ineq];
    X{i}.eq = [Xall.eq; X{i}.eq; W.eq];
end

%time range of states
if ~isfield(options.dynamics, 'Tmin') || isempty(options.dynamics.Tmin)
    options.dynamics.Tmin = zeros(nsys, 1);
    options.dynamics.Tmax = ones(nsys, 1)*options.Tmax;
end

if options.Tmax == Inf
    TIME_INDEP = 1;       
    t = [];
else    
    TIME_INDEP = 0;    
    assert(max(options.dynamics.Tmax) < Inf);
    if isempty(options.var.t)        
        t = sdpvar(1,1);
        options.var.t = t;
    else
        t = options.var.t;
    end
    nvar = nvar + 1;
end

%number of objectives (1 standard, 2+ minimum)
nobj = length(options.obj);




%% Form SOS constraints
%all variables
vars = [t; x; w];

%Set up the test function and dynamics
[v, cv] = polynomial(vars, d);
gamma = sdpvar(1, 1);


%accumulate all coefficients and constraints
cons = [];
coeff_list = [gamma; cv];

%initial condition Psatz
%v - gamma >= 0 on X0
if TIME_INDEP
    v0 = v;
else
    v0 = replace(v, t, 0);
end
vars0 = [x; w];

[p0, cons0, coeff0] = constraint_psatz(v0 - gamma, X0, vars0, d);

cons = [cons; cons0];
coeff_list = [coeff_list; coeff0];

%switching dynamics
for i = 1:nsys
    
    Lv = jacobian(v, x) * f{i};
    Xf = X{i};
    
    if ~TIME_INDEP
        Lv = Lv + jacobian(v, t);
        
        time_bound = (options.dynamics.Tmax(i) - t ) * (t - options.dynamics.Tmin(i));
        Xf.ineq = [Xf.ineq; time_bound];
    end
    [pf, consf, coefff] = constraint_psatz(Lv, Xf, vars, d);
    
    cons = [cons; consf];
    coeff_list = [coeff_list; coefff];
end


%deal with objective
%add scaling later
obj = options.obj;
if nobj == 1
    %single objective
    [pc, consc, coeffc] = constraint_psatz(-v - obj, Xall, vars, d);
    
    cons = [cons; consc];
    coeff_list = [coeff_list; coeffc];   
else
    %minimum of multiple objectives
    %beta is a simplex weighting between different objectives
    beta = sdpvar(nobj, 1);
    [pc, consc, coeffc] = constraint_psatz(-v - sum(beta.*obj), Xall, vars, d);
    
    cons = [cons; consc; sum(beta) == 1; beta >= 0];
    coeff_list = [coeff_list; coeffc; beta];
end

%objective of the dual problem. Yes this is confusing
if isempty(options.prev_cost)
    objective = -gamma;
else
    %feasibility program based on https://yalmip.github.io/Strictly-feasible-sum-of-squares/
    %but it doesn't work
    objective = [];
    cons = [cons; gamma >= options.prev_cost - 0.0001];
end

%% Solve program and extract results
%set up problem
opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;

%doesn't really help
% opts.sos.numblk =  1e-6;


[sol, monom, Gram, residual] = solvesos(cons, objective, opts, coeff_list);

peak_val = value(gamma);



%Yalmip reformulates the problem, therefore unable to extract moments from
%dual solution. This is unfortunate, but is the price for higher accuracy.

v_rec = value(cv)' * monolist(vars, d);

% %no more scaling should be necessary from this point on
% v = dual_rec'*monp_unscale;
Lv_rec = [];
for i = 1:nsys
    Lv_curr = jacobian(v_rec, x)*options.dynamics.f{i};
    if ~TIME_INDEP
        Lv_curr = Lv_curr + jacobian(v_rec, t);
    end
    Lv_rec = [Lv_rec; Lv_curr];
end

%% Output results to data structure
%out = 1;

out = struct;

%recover optima
out.order = order;
out.peak_val = peak_val;
out.monom = monom;
out.Gram = Gram;
out.mineig = min(cellfun(@(q) min(eig(q)), Gram));
out.feasible = out.mineig > 0;
out.residual = residual;
out.optimal = 0;

if nobj == 1
    out.beta = 1;
else
    out.beta = value(beta);
end


out.var = struct('t', t, 'x', x, 'w', w);


%Functions for simulating system dynamics

out.func = struct;
%polynomials
out.func.v = v_rec;
out.func.Lv = Lv_rec;

out.func.fval = cell(nsys, 1);  %dynamics
out.func.Xval = cell(nsys, 1);  %support set
out.func.event = cell(nsys, 1); %Modification for ode45 event
% 
for i = 1:nsys
    
    if nw > 0
        %fval_curr = @(t,x,w) replace(options.dynamics.f{i}, [tp; x; wp], [t; x; wp]);
        
        if TIME_INDEP
            fval_curr = @(tt,xt, wt) replace(options.dynamics.f{i}, [x; w], [xt; wt]);
        else
            fval_curr = @(tt,xt, wt) replace(options.dynamics.f{i}, [t; x; w], [tt; xt; wt]);
        end
        
%         
        Xval_curr = @(xt,wt) all(replace([options.dynamics.X{i}; XR], [x, w], [xt,wt]));
%         
%         %space is inside support, time between Tmin and Tmax
%         %if event_curr=0 stop integration, and approach from either
%         %direction (only going to be in negative direction, given that
%         %event is a 0/1 indicator function
%         event_curr = @(t,x,w) deal(all([Xval_curr(x,w); ...
%             t >= options.dynamics.Tmin(i); t < options.dynamics.Tmax(i)]), 1, 0);
    else       
        if TIME_INDEP
            fval_curr = @(tt,xt) replace(options.dynamics.f{i}, x, xt);
        else
            fval_curr = @(tt,xt) replace(options.dynamics.f{i}, [t; x], [tt; xt]);
        end
        Xval_curr = @(xt) constraint_eval(X{i}, x, xt);
        %Xval_curr = @(xt) all(replace([options.dynamics.X{i}; XR], x, xt));
        %Xval_curr = @(xt) all([replace([options.dynamics.X{i}.ineq; XR], x, xt) >=0;...
        %    abs(replace([options.dynamics.X{i}.eq], x, xt)) <=1e-8]);
        
        %event_curr = @(t,x) deal(all([Xval_curr(x); ...
        %    t >= options.dynamics.Tmin(i); t < options.dynamics.Tmax(i)]), 1, 0);
        event_curr = @(tt, xt) support_event(tt, xt, Xval_curr, ...
            options.dynamics.Tmin(i), options.dynamics.Tmax(i));
    end
%     
    out.func.fval{i} = fval_curr;
    out.func.Xval{i} = Xval_curr;        
    out.func.event{i} = event_curr;
end
% 

out.dynamics = struct;
out.dynamics.f = out.func.fval;
out.dynamics.event = out.func.event;

out.dynamics.time_indep = TIME_INDEP;


%functions that should be nonnegative along valid trajectories
if nobj > 1
    out.func.cost = @(xt) min(polyval_yalmip(options.obj, x, xt), 1);
else
    out.func.cost = @(xt) (polyval_yalmip(options.obj, x, xt));
end
 
% %TODO: missing the 'beta' weights for cost function in this representation
% %under multiple cost functions. Check that out later, proper dual
% %representation and verification of nonnegativity
if TIME_INDEP
    
%     out.func.vval = @(xt) replace(v_rec, x, xt);    %dual v(t,x,w)
%     out.func.Lvval = @(xt) replace(Lv_rec, x, xt);   %Lie derivative Lv(t,x,w)

    out.func.vval  = @(xt) polyval_yalmip(v_rec, x, xt);    %dual v(t,x,w)
    out.func.Lvval = @(xt) polyval_yalmip(Lv_rec, x, xt);   %Lie derivative Lv(t,x,w)


    %out.func.nonneg = @(xt) [out.func.vval(xt) - peak_val; out.func.Lvval(xt); -out.func.vval(xt) - out.func.cost(xt)];
    out.func.nonneg = @(xt) [out.func.vval(xt) - peak_val, out.func.Lvval(xt), -out.func.vval(xt) - out.func.cost(xt)]';
else
    out.func.vval  = @(tt, xt) polyval_yalmip(v_rec, [t; x], [tt; xt]);    %dual v(t,x,w)
    out.func.Lvval = @(tt, xt) polyval_yalmip(Lv_rec, [t; x], [tt; xt]);   %Lie derivative Lv(t,x,w)

%     out.func.nonneg = @(tt, xt) [out.func.vval(tt, xt) - peak_val; out.func.Lvval(t, x); -out.func.vval(tt, xt) - out.func.cost(xt)];
out.func.nonneg = @(tt, xt) [out.func.vval(tt, xt) - peak_val, out.func.Lvval(t, x), -out.func.vval(tt, xt) - out.func.cost(xt)]';
end
% 
% 
out.dynamics.nonneg = out.func.nonneg;
%% Done!

end

function pval = polyval_yalmip(poly, vars, pt)
%evaluate an sdpvar polynomial 'poly' in vars 'vars' at points 'pt'

npts = size(pt, 2);

pval = zeros(npts, 1);
for i = 1:npts
    ptcurr = pt(:, i);
    pval(i) = replace(poly, vars, ptcurr);
end

end