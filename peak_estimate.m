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
    order = 2;
end
d = 2*order;

%% Process option data structure

%state variable
xs = p_opt.var.x;
nx = length(xs);

%time variable
if ~isempty(p_opt.var.t)
    ts = p_opt.var.t;
else
    ts = 0;
end

%parameter variables
if ~isempty(p_opt.var.w)
    PARAM = 1;
    ws = p_opt.var.w;
    %Wsupp_f = matlabFunction(p_opt.param);
else
    PARAM = 0;
end

%Time independent measures if valid
if isempty(p_opt.var.t) || p_opt.Tmax = Inf
    TIME_INDEP = 1;
else
    TIME_INDEP = 0;
end


%number of switching subsystems
if iscell(p_opt.dynamics.f)
    nsys = length(p_opt.dynamics.f);
else
    nsys = 1;
end

%number of objectives (1 standard, 2+ minimum)
nobj = length(p_opt.obj);

%number of times at which the space support is specified
if iscell(p_opt.state_fix.X)
    nbreaks = length(p_opt.state_fix.X);
    time_breaks = p_opt.state_fix.T;
else
    nbreaks = 1;
end


%% Measures and variables
%measures in the problem:
%start with no time breaks, one initial measure
%one occupation measure per switched systems
%one peak measure
x_occ = mpol('xocc', [nx, nsys]);

xp = mpol('xp', [nx, 1]);

%replace with time breaks
x0 = mpol('x0', [nx, 1]);

mu = cell(nsys, 1);
v = cell(nsys, 1);

if TIME_INDEP
           
    mup = meas(xp);
    mu0 = meas(x0);
    [mu, X, Ay] = occupation_measure(f, X, var, x_occ, d)
    for i = 1:nsys
        xcurr = x_occ(:, i);
        mu{i} = meas(xcurr);
        v{i}  = mmon(xcurr, d);
    end
else
    %not actually sure, maybe all switching occupation measures have the same titme
    t_occ = mpol('tocc', [1, nsys]);
    tp = mpol('tp', 1);
    t0 = mpol('t0', 1);
    
    mup = meas([tp; xp]);
    mu0 = meas([t0; x0]);
    
    for i = 1:nsys
        xcurr = x_occ(:, i);
        mu{i} = meas([t_occ(i), xcurr]);
        v{i}  = mmon([t_occ(i), xcurr], d);
    end
    
end





%% Solve program and extract results


%% Output results to data structure
out = 1;

% out = struct;
% out.peak_val = peak_val;
% out.optimal = optimal;
% out.x0 = x0_rec;
% out.xp = xp_rec;
% out.Mfix = {};
% out.Mp = [];
% 
% out.v = 0;
% out.Lv = 0;

%% Done!

end

function [mu, X_occ, Ay] = occupation_measure(f, X, var, x_occ, d)
%form the occupation measure
%Input:
%   f:      Dynamics (polynomial)
%   X:      Support set upon which dynamics take place
%   var:    variables of problem
%   x_occ:  New variables for occupation measure
%   d:      Degree of measure
%
%Output:
%   mu:     Measure
%   X:      Support of measure
%   Ay:     Adjoint of lie derivative, Liouville

    var_all = [var.t; x_occ; var.w];

    f_occ = subs(f, var.x, x_occ);
    
    mu = meas(var_all);
    v = mmon(var_all, d);

    Ay = mom(diff(v, x_occ)*f_occ);
    if ~isempty(var.t)
        Ay = Ay + mom(diff(v, var.t));
    end
    
    X_occ = subs(X, var.x, x_occ);

end

