function [out] = peak_estimate(options, order)
%PEAK_ESTIMATE  Find the maximum value of a function p(x) along all
%trajectories starting from an initial set X0 along polynomial dynamics
%f(t, x). This is for finite time [0, Tmax], but if f(t, x) is autonomous
%then the problem can be solved for infinite time.
%   Accomodates switching as well.
%
%Input: options structure (peak_options) with fields:
%   var:        Structure of symbolic variables
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
%   state_fix:  Supports of trajectories at specific time instances
%       List:       Initial Support X0
%       Struct:     Fields X and T. X: support of trajectories at times T
%       (will need a better name. Generalizes X0)
%               X is a cell of lists of polynomial. 
%               T is an array ofdoubles             
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
%   xp:         Optimal point (obj(xp) = peak_val)
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
mset clearmeas

%state variable
nx = length(options.var.x);
nvar = nx;

%time variable
if ~isempty(options.var.t)
    nvar = nvar + 1;
    %TIME_INDEP = 0;
else
    %TIME_INDEP = 1;
end

%parameter variables
nw = length(options.var.w);
nvar = nvar + nw;

%compact support (artificial)
XR = (sum(options.var.x.^2) <= options.R^2);

%number of switching subsystems
if iscell(options.dynamics.f)
    nsys = length(options.dynamics.f);
    for i = 1:nsys
        options.dynamics.X{i} = [options.dynamics.X{i}; XR];
    end
else
    nsys = 1;
    options.dynamics.f = {options.dynamics.f};
    options.dynamics.X = {[options.dynamics.X; XR]};    
end

%time range of states
if ~isfield(options.dynamics, 'Tmin') || isempty(options.dynamics.Tmin)
    options.dynamics.Tmin = zeros(nsys, 1);
    options.dynamics.Tmax = ones(nsys, 1)*options.Tmax;
end

if options.Tmax == Inf
    TIME_INDEP = 1;       
else    
    TIME_INDEP = 0;    
    
    if isempty(options.var.t)
        mpol('t', 1, 1)
        tp = t;
        options.var.t = tp;
    else
        tp = options.var.t;
    end
    nvar = nvar + 1;
end


%number of objectives (1 standard, 2+ minimum)
nobj = length(options.obj);

%number of times at which the space support is specified
%do this later
% if iscell(options.state_fix.X)
%     nbreaks = length(options.state_fix.X);
%     time_breaks = options.state_fix.T;
% else
%     nbreaks = 1;
% end


%% Measures and variables
%measures in the problem:
%start with no time breaks, one initial measure
%one occupation measure per switched systems

%one peak measure
xp = options.var.x;
Xp = [options.state_supp; XR];

%deal with hanging variables and measures by letting the original x be the
%peak measure

%replace with time breaks
mpol('x0', nx, 1);
X0 = subs([options.state_init; XR], options.var.x, x0);



if nw > 0
    %parameters
    wp = options.var.w;    
    mpol('w0', nw, 1);
    
    Wp = options.param;
    W0 = subs(Wp, wp, w0);
    
    W = [Wp; W0; wp == w0];
    mpol('w_occ', nw, nsys);   
else
    W = [];    
    w_occ = zeros(0, nsys);
    wp = [];
    w0 = [];
end

%occupation measure
mpol('x_occ', nx, nsys);



mu = cell(nsys, 1);
v = cell(nsys, 1);
Ay = 0;
%X = cell(nsys, 1);
X_occ = []; %support

%measure information
if TIME_INDEP           
    mup = meas([xp; wp]);
    tp = [];
    monp = mmon([xp; wp], d);
    yp = mom(monp);
    
    mu0 = meas([x0; w0]);
    t0 = [];
    mon0 = mmon([x0; w0], d);
    y0 = mom(mon0);
    
    %[mu, X, Ay] = occupation_measure(f, X, options.var, x_occ, d)
    
    for i = 1:nsys
        %xcurr = x_occ(:, i);
        var_new = struct('t', [], 'x', x_occ(:, i), 'w', w_occ(:, i));
        
        [mu{i}, X_occ_curr, Ay_curr] = occupation_measure(options.dynamics.f{i}, ...
            options.dynamics.X{i}, options.var, var_new, d);
        
        X_occ = [X_occ; X_occ_curr];
        Ay = Ay + Ay_curr;
        
        if nw > 0
            W_curr = subs(Wp, wp, w_occ(:, i));
            W = [W; W_curr];
        end
    end
else
    %define time variables
    mpol('t_occ', 1, nsys);        
    
    %peak time
    %mpol('tp', 1);
    Xp = [Xp; tp*(options.Tmax - tp) >= 0];
    mup = meas([tp; xp; wp]);
    monp = mmon([tp; xp; wp], d);
    yp = mom(monp);
    
    %initial time
    mpol('t0', 1);
    %assign(t0, 0);
    X0 = [X0; t0 == 0];                
    mu0 = meas([t0; x0; w0]); 
    mon0 = mmon([t0; x0; w0], d);
    y0 = mom(mon0);
    
    
    for i = 1:nsys
        var_new = struct('t', t_occ(i), 'x', x_occ(:, i), 'w', w_occ(:, i));
        xcurr = x_occ(:, i);
        mu{i} = meas([t_occ(i); xcurr]);
        mon{i}  = mmon([t_occ(i); xcurr], d);
        
        [mu{i}, X_occ_curr, Ay_curr] = occupation_measure(options.dynamics.f{i}, ...
            options.dynamics.X{i}, options.var, var_new, d);
        
        %support the valid time range
        Tmin_curr = options.dynamics.Tmin(i);
        Tmax_curr = options.dynamics.Tmax(i);
        
        t_cons = (t_occ(i) - Tmin_curr)*(Tmax_curr - t_occ(i));            
        
        X_occ = [X_occ; X_occ_curr; t_cons >= 0];
        Ay = Ay + Ay_curr;
        
        if nw > 0
            W_curr = subs(Wp, wp, w_occ(:, i));
            W = [W; W_curr];
            %is 'w_occ(:, i) == wp' enforced implicitly?
        end
    end
    
end

%% Form Constraints and Solve Problem

supp_con = [X0; Xp; X_occ; W];
%careful of monic substitutions ruining dual variables
mom_con = [mass(mu0) == 1; Ay + y0 == yp];

if nobj == 1
    cost = subs(options.obj, [options.var.t; options.var.x], [tp; xp]);
else
    %add a new measure for the cost
    %honestly I only want a free variable, so bounds on support are not
    %necessary. The free variable is the off-diagonal entry of a 2x2 psd
    %matrix with top corner 1. Inefficient, but that's how gloptipoly is
    %interfaced.
    
    mpol('c');
    muc = meas(c);
    momc = mom(c);
    cost = momc;
    
    for i = 1:nobj
        curr_obj = subs(options.obj(i), var.x, xp);
        curr_mom_con = (momc <= mom(curr_obj));
        mom_con = [mom_con; curr_mom_con];        
    end       
end

objective = max(cost);

%% Solve program and extract results
%set up problem
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));

P = msdp(objective, ...
    mom_con, supp_con);



%solve LMI moment problem    

%for the flow problem, peak_test_alt has 1 substitution
%this script has 2 moment substitutions, so dual_rec does not have the same
%size as the vernoese map vv. I can't see a difference. What is going on?
[status,obj_rec, m,dual_rec]= msol(P);

M0 = double(mmat(mu0));
Mp = double(mmat(mup));
if nobj > 1
    Mc = double(mmat(muc));
end
%1, t, x1, x2
M0_1 = M0(1:(nvar+1), 1:(nvar+1));
Mp_1 = Mp(1:(nvar+1), 1:(nvar+1));

rank0 = rank(Mp_1, options.rank_tol);
rankp = rank(M0_1, options.rank_tol);



x0_rec = double(mom(x0));
xp_rec = double(mom(xp));
% 
% if ~TIME_INDEP    
%     tp_out = T*double(mom(tp));
% end

%set up dual variable
%syms xc [nx, 1]
%scaling with t



%vv = monolist([tc/T; xc; yc], d);
%vv = conj(monolistYalToGlop([tc/T; xc; yc], d));
% if TIME_INDEP
%     vv = conj(monolistYalToGlop(xc, d));
% else
%     syms tc
%     vv = conj(monolistYalToGlop([tc; xc], d));    
% end

%recovered dual variables from msdp, correspond to the free
%variables in the sedumi problem

%set up polynomials
% p = dual_rec'*vv;
% Lp = cell(nsys, 1);
% for i = 1:nsys
%     Lp{i} = jacobian(p, xc)*options.dynamics.f{i}(xc);
% end
% 
% %turn into functions
% pval = matlabFunction(p);
% if nsys == 1
%     Lpval = matlabFunction(Lp{1});
% else
%     Lpval = cellfun(@(Lpf) matlabFunction(Lpf), Lp, 'UniformOutput', false);
% end
v = dual_rec'*monp;
Lv = [];
for i = 1:nsys
    Lv_curr = diff(v, xp)*options.dynamics.f{i};
    if ~TIME_INDEP
        Lv_curr = Lv_curr + diff(v, tp);
    end
    Lv = [Lv; Lv_curr];
end

%% Output results to data structure
%out = 1;

out = struct;

%recover optima
out.peak_val = -obj_rec;
out.optimal = (rank0 == 1) && (rankp == 1);
out.x0 = x0_rec;
out.xp = xp_rec;
out.M0 = M0_1;
out.Mp = Mp_1;

if ~TIME_INDEP
    out.tp = double(mom(tp));
end

out.var = struct('t', tp, 'x', xp, 'w', wp);

%functions and dual variables
out.func = struct;
out.func.dual_rec = dual_rec;
out.func.v = v;
out.func.Lv = Lv;

%functions that should be nonnegative along valid trajectories
out.func.cost = @(x) min(eval(options.obj, xp, x));
out.func.vval = @(x) eval(v, xp, x);    %dual v(t,x,w)
out.func.Lvval = @(x) eval(Lv, xp, x);   %Lie derivative Lv(t,x,w)

%TODO: missing the 'beta' weights for cost function in this representation
%under multiple cost functions. Check that out later, proper dual
%representation and verification of nonnegativity
out.func.nonneg = @(x) [out.func.vval(x) + obj_rec; out.func.Lvval(x); -out.func.vval(x) - out.func.cost(x)];

%evaluate system points
out.func.fval = cell(nsys, 1);  %dynamics
out.func.Xval = cell(nsys, 1);  %support set
out.func.event = cell(nsys, 1); %Modification for ode45 event

for i = 1:nsys
    
    if nw > 0
        %fval_curr = @(t,x,w) eval(options.dynamics.f{i}, [tp; xp; wp], [t; x; wp]);
        
        if TIME_INDEP
            fval_curr = @(t,x, w) eval(options.dynamics.f{i}, [xp; wp], [x; w]);
        else
            fval_curr = @(t,x, w) eval(options.dynamics.f{i}, [tp; xp; wp], [t; x; w]);
        end
        
        
        Xval_curr = @(x,w) all(eval(options.dynamics.X{i}, [xp, w], [x,w]));
        
        %space is inside support, time between Tmin and Tmax
        %if event_curr=0 stop integration, and approach from either
        %direction (only going to be in negative direction, given that
        %event is a 0/1 indicator function
        event_curr = @(t,x,w) deal(all([Xval_curr(x,w); ...
            t >= options.dynamics.Tmin(i); t < options.dynamics.Tmax(i)]), 1, 0);
    else       
        if TIME_INDEP
            fval_curr = @(t,x) eval(options.dynamics.f{i}, xp, x);
        else
            fval_curr = @(t,x) eval(options.dynamics.f{i}, [tp; xp], [t; x]);
        end
        Xval_curr = @(x) all(eval(options.dynamics.X{i}, xp, x));
        
        %event_curr = @(t,x) deal(all([Xval_curr(x); ...
        %    t >= options.dynamics.Tmin(i); t < options.dynamics.Tmax(i)]), 1, 0);
        event_curr = @(t, x) support_event(t, x, Xval_curr, ...
            options.dynamics.Tmin(i), options.dynamics.Tmax(i));
    end
    
    out.func.fval{i} = fval_curr;
    out.func.Xval{i} = Xval_curr;        
    out.func.event{i} = event_curr;
end

out.dynamics = struct;
out.dynamics.f = out.func.fval;
out.dynamics.event = out.func.event;


%% Done!

end

function [mu, X_occ, Ay] = occupation_measure(f, X, var, var_new, d)
%form the occupation measure
%Input:
%   f:      Dynamics (function)
%   X:      Support set upon which dynamics take place
%   var:    variables of problem
%   x_occ:  New variables for occupation measure
%   d:      Degree of measure
%
%Output:
%   mu:     Measure
%   X:      Support of measure
%   Ay:     Adjoint of lie derivative, Liouville

    %var_all = [var.t; x_occ; var.w];

    x_occ = var_new.x;
    
    if ~isempty(var_new.t)
        t_occ = var_new.t;
    else
        t_occ = [];
    end
    
    if ~isempty(var_new.w)
        w_occ = var_new.w;
    else
        w_occ = [];
    end
    
    vars_prev = [var.t; var.x; var.w];
    vars_new = [t_occ; x_occ; w_occ];
    
    %f_occ = subs(f, var.x, x_occ);
    f_occ = subs(f, vars_prev, vars_new);
    %f_occ = f(var_all);
    
    mu = meas(vars_new);
    mon = mmon(vars_new, d);

    Ay = mom(diff(mon, x_occ)*f_occ);
    if ~isempty(var.t)
        Ay = Ay + mom(diff(mon, t_occ));
    end
    
    X_occ = subs(X, var.x, x_occ);

end

