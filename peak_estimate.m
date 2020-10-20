function [out] = peak_estimate(options, order)
%PEAK_ESTIMATE  Find the maximum value of a function p(x) along all
%trajectories starting from an initial set X0 along polynomial dynamics
%f(t, x). This is for finite time [0, Tmax], but if f(t, x) is autonomous
%then the problem can be solved for infinite time.
%   Accomodates switching as well.
%
%Input: options structure (peak_options) with fields:
%   var:        Structure of symbolic variables (@mpol)
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
%   state_supp: Support set of total set X (@supcon)
%   state_init: Support set of initial set X0 (@supcon)
%   param:      Support of parametric uncertainty w (if var.w nonempty)
%       
%   rank_tol:   Rank tolerance for moment matrix to be rank-1
%   box:        Box containing valid region X, default to [-1, 1]^n
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
%   M0:         Moment Matrix of Initial Distribution
%   Mp:         Moment Matrix of Peak Distribution
%   Mocc:       Moment Matrix of Occupation Measure
%
%   v:          Polynomial from dual program, v(t, x) - gamma >=0 invariant
%   Lv:         Lie derivatives along trajectories of system for v
%   beta:       multiple cost tradeoff at optimal point

if nargin < 2
    order = 3;
end
degree = 2*order;

%% Process option data structure
mset clearmeas

%state variable
nx = length(options.var.x);
nvar = nx;

%parameter variables
nw = length(options.var.w);
nvar = nvar + nw;


nd = length(options.var.d);
nvar = nvar + nd;

nb = length(options.var.b);
nvar = nvar + nb;

%compact support (artificial)
%XR = (sum(options.var.x.^2) <= options.R^2);
%scaling of the box
xp = options.var.x;
%TODO: deal with parameters nx -> nx + nw;
if isempty(options.box) || (length(options.box)==1 && options.box == 0)
    XR_scale = [];
    XR_unscale = [];
    
    xp_scale = xp;
    xp_inv_scale = xp;
else    
    [box, box_center, box_half] = box_process(nx, options.box);

    XR_scale = (xp.^2) <= 1;
    XR_unscale = (xp - box_center).^2 <= box_half.^2;
    
    xp_scale = box_half.*xp + box_center;
    xp_inv_scale = (xp - box_center) .* (1./box_half);
end
Xsupp = options.state_supp;

%XR_scale = [xp <= 1; xp >= -1];
%XR_unscale = [xp - box_center <= box_half; xp - box_center >= -box_half];


if options.scale
    XR = XR_scale;
else
    XR = XR_unscale;
end



%number of switching subsystems
f = options.dynamics.f;
X = options.dynamics.X;

if ~isfield(options.dynamics, 'discrete')
    options.dynamics.discrete = 0;
end

if iscell(options.dynamics.f)
    nsys = length(options.dynamics.f);
else
    nsys = 1;
    options.dynamics.f = {options.dynamics.f};
    options.dynamics.X = {options.dynamics.X};    
    f = {f};
    X = {X};
end

%time range of dynamics
if ~isfield(options.dynamics, 'Tmin') || isempty(options.dynamics.Tmin)
    options.dynamics.Tmin = zeros(nsys, 1);
    options.dynamics.Tmax = ones(nsys, 1)*options.Tmax;
end


%time range of all trajectories
if options.Tmax == Inf || options.dynamics.discrete
    TIME_INDEP = 1;    
    Tsupp = [];
    tp = [];
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
    
    Tsupp = tp*(1-tp) >= 0;

end



%now scale everything
if options.scale    
    for i = 1:nsys
        f{i} = subs((1./box_half) .*f {i}, xp, xp_scale);
        
        X{i} = [subs([X{i}; Xsupp], xp, xp_scale); XR_scale];
    end
else
    for i = 1:nsys
       X{i} = [X{i}; Xsupp; XR_unscale];
    end
end

if ~TIME_INDEP
    for i = 1:nsys
        f{i} = options.Tmax * subs(f{i}, tp, tp*options.Tmax);
    end
end

%number of objectives (1 standard, 2+ minimum)
nobj = length(options.obj);

%% Measures and variables
%measures in the problem:
%start with no time breaks, one initial measure
%one occupation measure per switched systems

%one peak measure

if options.scale
%     Xp = [subs(options.state_supp, xp, xp_scale); Xsupp; XR];
    Xp = [Tsupp; subs(options.state_supp, xp, xp_scale); XR]; %that might be the bug (hopefully)
else
    
    %this is where dual_rec breaks
%     Xp = [Xsupp; XR];
    %adding equality constraints here breaks the invariant function
    %dual_rec'monomials. Since X0 and X_occ are constrained within Xsupp,
    %this should not be a problem. If dynamics stay within Xsupp, then Xp
    %should be supported in Xsupp at optimality.
    Xp = [Tsupp; Xsupp; XR];
end
%deal with hanging variables and measures by letting the original x be the
%peak measure

%replace with time breaks
mpol('x0', nx, 1);

if options.scale
    X0_scale = subs(options.state_init, xp, xp_scale);
    
    X0 = subs([X0_scale; Xsupp; XR], xp, x0);
else
    X0 = subs([options.state_init; Xsupp; XR], xp, x0);
end


%deal with time-independent uncertainty (parameters)
if nw > 0
    %parameters
    wp = options.var.w;    
    mpol('w0', nw, 1);
    
    Wp = options.param;
%     W0 = subs(Wp, wp, w0);
    W = Wp;
    %Wp = W0 since w is time-independent
%     W = [Wp; W0];  
    
else
    W = [];    
    w_occ = zeros(0, nsys);
    wp = [];
    w0 = [];
end

if TIME_INDEP
    t0 = [];
else
    mpol('t0', 1, 1);
    X0 = [t0 == 0; X0];
end


mu_occ = cell(nsys, 1);
mu_occ_sum = 0;
mon_occ = cell(nsys, 1);
v = cell(nsys, 1);
Ay = 0;
%X = cell(nsys, 1);
X_occ = []; %support

%% occupation measure
%state
mpol('x_occ', nx, nsys);

%time-independent parameter
if nw > 0
    mpol('w_occ', nw, nsys);
else
    w_occ = zeros(0, nsys);
end

%time-dependent general uncertainty
if nd > 0
    mpol('d_occ', nd, nsys);
else
    d_occ = zeros(0, nsys);
end
D = [];

%time-dependent box in [0, 1]
if nb > 0
    mpol('b_occ', nb, nsys);
else
    b_occ = zeros(0, nsys);
end
%this will be tricky. deal with it later.
%involves box control-occupation measures and absolute continuity
%constraints.


%measure information
if TIME_INDEP           
    mup = meas([xp; wp]);
    tp = [];
    monp = mmon([xp; wp], degree);
    yp = mom(monp);
    
    mu0 = meas([x0; w0]);
    t0 = [];
    mon0 = mmon([x0; w0], degree);
    y0 = mom(mon0);        
    
    for i = 1:nsys
        var_new = struct('t', [], 'x', x_occ(:, i), 'w', w_occ(:, i), ...
            'd', d_occ(:, i), 'b', b_occ(:, i));
        
        [mu_occ{i}, mon_occ{i}, X_occ_curr, Ay_curr] = occupation_measure(f{i}, ...
            X{i}, options.var, var_new, degree, options.dynamics.discrete);
        
        %mu_occ_sum = mu_occ_sum + mu_occ{i};
        X_occ = [X_occ; X_occ_curr];
        Ay = Ay + Ay_curr;
        
        if nw > 0
            W_curr = subs(Wp, wp, w_occ(:, i));
            W = [W; W_curr];
        end
        
        if nd > 0
            D_curr = subs(options.disturb, options.var.d, d_occ(:, i));
            D = [D; D_curr];
        end
    end
else
    %define time variables
    mpol('t_occ', 1, nsys);                    
    
    %peak time
    %mpol('tp', 1);
    mup = meas([tp; xp; wp]);
    monp = mmon([tp; xp; wp], degree);
    yp = mom(monp);
    
    %initial time             
    mu0 = meas([t0; x0; w0]); 
    mon0 = mmon([t0; x0; w0], degree);
    y0 = mom(mon0);
    
    
    for i = 1:nsys
        var_new = struct('t', t_occ(i), 'x', x_occ(:, i), 'w', w_occ(:, i),...
            'd', d_occ(:, i), 'b', b_occ(:, i));
        
        [mu_occ{i}, mon_occ{i}, X_occ_curr, Ay_curr] = occupation_measure(...
            f{i}, X{i}, options.var, var_new, degree, 0);
        
        %support the valid time range
        Tmin_curr = options.dynamics.Tmin(i);
        Tmax_curr = options.dynamics.Tmax(i);
        
        t_cons = (t_occ(i) - Tmin_curr/options.Tmax)*(Tmax_curr/options.Tmax - t_occ(i));            

        %mu_occ_sum = mu_occ_sum + mu_occ{i};
        X_occ = [X_occ; t_cons >= 0; X_occ_curr];
        Ay = Ay + Ay_curr;
        
        if nw > 0
            W_curr = subs(Wp, wp, w_occ(:, i));
            W = [W; W_curr];            
        end
        
        if nd > 0
            D_curr = subs(options.disturb, options.var.d, d_occ(:, i));
            D = [D; D_curr];
        end
    end
    
end

%% Form Constraints and Solve Problem

%supp_con = [Xp; X0; X_occ; W];
supp_con = [Xp; X0; X_occ; W; D];
%careful of monic substitutions ruining dual variables
Liou = Ay + (y0 - yp);
mom_con = [mass(mu0) == 1; Liou == 0];

obj = options.obj;
if nobj == 1
    if TIME_INDEP
        cost = obj;
    else
        cost = subs(obj, tp, tp*options.Tmax);
    end
    if options.scale
        cost = subs(obj, xp, xp_scale);
    end
    %cost = subs(options.obj, [options.var.t; options.var.x], [tp; xp]);
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
    
%     cost_mom_con = (mass(muc) == 1);
    cost_mom_con = (mass(muc) == 1);
%     cost_mom_con = [];
    for i = 1:nobj
        %curr_obj = subs(options.obj(i), var.x, xp);
        if TIME_INDEP   
            curr_obj = options.obj(i);
        else
            curr_obj = subs(options.obj(i), tp, tp*options.Tmax);
        end
        if options.scale
            curr_obj = subs(options.obj(i), xp, xp_scale);
        end
        
        cost_mom_con = [cost_mom_con; (momc <= mom(curr_obj))];
        
    end 
    
    mom_con = [cost_mom_con; mom_con];
end



if isempty(options.prev_cost)
    objective = max(cost);
else
    objective = 0;
    %should this be a support or a moment?
    supp_con = [supp_con; cost >= 0.9999*options.prev_cost];
end

%% Solve program and extract results
%set up problem
mset('yalmip',true);
mset(sdpsettings('solver', options.solver));

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

rank0 = rank(M0_1, options.rank_tol);
rankp = rank(Mp_1, options.rank_tol);

x0_rec = double(mom(x0));
xp_rec = double(mom(xp));

if ~TIME_INDEP
    tp_rec = double(mom(tp));
end

%how does this work with moment substitutions?
%necessary with c^2 + s^2 = 1
if options.scale
    x0_rec = eval(xp_inv_scale, xp, x0_rec);
    xp_rec = eval(xp_inv_scale, xp, xp_rec);
    monp_unscale = subs(monp, xp, xp_inv_scale);
else
    monp_unscale = monp;
end

if ~TIME_INDEP  
    tp_rec = tp_rec * options.Tmax;
    monp_unscale = subs(monp, tp, tp/options.Tmax);        
end

%no more scaling should be necessary from this point on

%with equality constraints, dual_rec is bigger than just v
%find the ordering of dual_rec, and which entries correspond to v
%v = dual_rec'*monp_unscale;
dual_rec_v = dual_rec{1};
% dual_rec_v = dual_rec_v((end - length(monp_unscale)+1):end);

dual_rec_v = dual_rec_v(1:length(monp_unscale));

if nobj == 1
    beta = 1;
else
    %this is an artifact of gloptipoly not having standard scalar
    %variables. A degree-1 measure is required.
    beta = dual_rec{2};
end

v = dual_rec_v'*monp_unscale;
Lv = [];
for i = 1:nsys
    if options.dynamics.discrete
        %discrete time
        pushforward = subs(v, xp, options.dynamics.f{i});
        Lv_curr = pushforward - v;
    else
        %continuous time
        Lv_curr = diff(v, xp)*options.dynamics.f{i};
        if ~TIME_INDEP
            Lv_curr = Lv_curr + diff(v, tp);
        end
    end
    Lv = [Lv; Lv_curr];
end

%% Output results to data structure
%out = 1;

out = struct;

%recover optima
out.order = order;
out.peak_val = obj_rec;
out.optimal = (rank0 == 1) && (rankp == 1);
out.x0 = x0_rec;
out.xp = xp_rec;
% out.M0 = M0_1;
% out.Mp = Mp_1;

%moment matrices
out.M0 = M0;
out.Mp = Mp;

%occupation measure
Mocc_cell=cellfun(@(m) double(mmat(m)), mu_occ, 'UniformOutput', false);
% Mocc = sum(cat(3, Mocc_cell{:}), 3);
% out.Mocc = Mocc;

out.Mocc = Mocc_cell;

%objective
if nobj > 1
    out.Mc = double(mmat(muc));
end

if ~TIME_INDEP
    out.tp = tp_rec;
end

if nw > 0
    out.w = double(mom(wp));
else
    out.w = [];
end

out.var = struct('t', tp, 'x', xp, 'w', wp);


%Functions for evaluating system dynamics
out.func.fval = cell(nsys, 1);  %dynamics
out.func.fval_all = cell(nsys, 1);  %dynamics
out.func.Xval = cell(nsys, 1);  %support set
out.func.event = cell(nsys, 1); %Modification for ode45 event

for i = 1:nsys
    Xval_curr = @(x) all(eval([options.dynamics.X{i}; XR_unscale], xp, x));
    
    
    
    %dynamics with all variables
    %useful for updated switch sampling code including time-varying
    %uncertainty
    if TIME_INDEP
        fval_curr_all = @(t, x, w, d, b) eval(options.dynamics.f{i}, ...
            [xp; wp; options.var.d; options.var.b], [x; w; d; b]);
    else
        fval_curr_all = @(t,x, w, d, b) eval(options.dynamics.f{i}, ...
            [tp; xp; wp; options.var.d; options.var.b], [t; x; w; d; b]);
    end
    
    if nw > 0
        %fval_curr = @(t,x,w) eval(options.dynamics.f{i}, [tp; xp; wp], [t; x; wp]);
        
        if TIME_INDEP
            fval_curr = @(t, x, w) eval(options.dynamics.f{i}, [xp; wp], [x; w]);
        else
            fval_curr = @(t,x, w) eval(options.dynamics.f{i}, [tp; xp; wp], [t; x; w]);
        end
        %space is inside support, time between Tmin and Tmax
        %if event_curr=0 stop integration, and approach from either
        %direction (only going to be in negative direction, given that
        %event is a 0/1 indicator function
        event_curr = @(t,x,w) support_event(t, x, Xval_curr, ...
            options.dynamics.Tmin(i), options.dynamics.Tmax(i)); %should w be present here?
    else       
        if TIME_INDEP
            fval_curr = @(t,x) eval(options.dynamics.f{i}, xp, x);
        else
            fval_curr = @(t,x) eval(options.dynamics.f{i}, [tp; xp], [t; x]);
        end
        event_curr = @(t, x) support_event(t, x, Xval_curr, ...
            options.dynamics.Tmin(i), options.dynamics.Tmax(i));
    end
    
    out.func.fval_all{i} = fval_curr_all;
    out.func.fval{i} = fval_curr;
    out.func.Xval{i} = Xval_curr;        
    out.func.event{i} = event_curr;
end

out.dynamics = struct;
out.dynamics.f = out.func.fval;
out.dynamics.f_all = out.func.fval_all;
out.dynamics.event = out.func.event;
out.dynamics.discrete = options.dynamics.discrete;
out.dynamics.time_indep = TIME_INDEP;
event_all = @(tt, xt) cell2mat(cellfun(@(ex) ex(tt, xt), out.func.event, 'UniformOutput', false));

%% functions and dual variables
out.func = struct;
out.func.dual_rec = dual_rec;
out.func.v = v;
out.func.beta = beta;
out.func.Lv = Lv;

%functions that should be nonnegative along valid trajectories
out.func.cost_all = cell(nobj, 1);
for i = 1:nobj
    out.func.cost_all{i} = @(x) (eval(options.obj(i), xp, x));
end

if nobj > 1
    out.func.cost = @(x) min(eval(options.obj, xp, x));
else
    out.func.cost = @(x) (eval(options.obj, xp, x));
end

%TODO: missing the 'beta' weights for cost function in this representation
%under multiple cost functions. Check that out later, proper dual
%representation and verification of nonnegativity

if nw > 0
    if TIME_INDEP    
        out.func.vval = @(x, w) eval(v, [xp; wp], [x; w*ones(1, size(x, 2))]);    %dual v(t,x,w)
        out.func.Lvval = @(x, w) eval(Lv, [xp; wp], [x; w*ones(1, size(x, 2))]);  %Lie derivative Lv(t,x,w)
        
        out.func.nonneg = @(x, w) [out.func.vval(x, w) + obj_rec; out.func.Lvval(x, w); -out.func.vval(x, w) - beta'*out.func.cost(x)];
%         out.func.nonneg = @(x, w) [out.func.vval(x, w) + obj_rec; out.func.Lvval(x, w).*event_all(zeros(1, size(x, 2)),x); -out.func.vval(x, w) - beta'*out.func.cost(x)];
    else
        out.func.vval = @(t, x, w) eval(v, [tp; xp; wp], [t; x; w*ones(1, size(x, 2))]);   %dual v(t,x,w)
        out.func.Lvval = @(t, x, w) eval(Lv, [tp; xp; wp],  [t; x; w*ones(1, size(x, 2))]);   %Lie derivative Lv(t,x,w)
        
        out.func.nonneg = @(t, x, w) [out.func.vval(t, x, w) + obj_rec; out.func.Lvval(t, x, w); -out.func.vval(t, x, w) - beta'*out.func.cost(x)];
%         out.func.nonneg = @(t, x, w) [out.func.vval(t, x, w) + obj_rec; out.func.Lvval(t, x, w).*event_all(t, x); -out.func.vval(t, x, w) - beta'*out.func.cost(x)];
    end

else
    if TIME_INDEP    
        out.func.vval = @(x) eval(v, xp, x);    %dual v(t,x,w)
        out.func.Lvval = @(x) eval(Lv, xp, x);   %Lie derivative Lv(t,x,w)

        out.func.nonneg = @(x) [out.func.vval(x) + obj_rec; out.func.Lvval(x); -out.func.vval(x) - out.func.cost(x)];
%         out.func.nonneg = @(x) [out.func.vval(x) + obj_rec; out.func.Lvval(x).*event_all(zeros(1, size(x, 2)),x); -out.func.vval(x) - out.func.cost(x)];
    else
        out.func.vval = @(t, x) eval(v, [tp; xp], [t; x]);    %dual v(t,x,w)
        out.func.Lvval = @(t, x) eval(Lv, [tp; xp], [t; x]);   %Lie derivative Lv(t,x,w)
        out.func.nonneg = @(t, x) [out.func.vval(t, x) + obj_rec; out.func.Lvval(t, x); -out.func.vval(t, x) - out.func.cost(x)];
%         out.func.nonneg = @(t, x) [out.func.vval(t, x) + obj_rec; out.func.Lvval(t, x).*event_all(t, x); -out.func.vval(t, x) - out.func.cost(x)];
    end
end

out.dynamics.cost = out.func.cost;
out.dynamics.nonneg = out.func.nonneg;
%% Done!

end

function [mu, mon, X_occ, Ay] = occupation_measure(f, X, var, var_new, degree, discrete)
%form the occupation measure
%Input:
%   f:        Dynamics (function)
%   X:        Support set upon which dynamics take place
%   var:      variables of problem [t,x,w,d,b]
%   var_new:  New variables for occupation measure
%   degree:   Degree of measure
%   discrete: True (discrete system) or False (continuous system)
%
%Output:
%   mu:     Occupation Measure
%   mon:    monomoials in Liouville constraint
%   X:      Support of measure in state
%   Ay:     Adjoint of lie derivative, Liouville

    %var_all = [var.t; x_occ; var.w];

    if nargin < 6
        discrete = 0;
    end
    
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
    
    if ~isempty(var_new.d)
        d_occ = var_new.d;
    else
        d_occ = [];
    end
    
    %box variables?
    
    vars_prev = [var.t; var.x; var.w; var.d];
    vars_new = [t_occ; x_occ; w_occ; d_occ];
    
    f_occ = subs(f, vars_prev, vars_new);
    
    mu = meas([t_occ; x_occ; w_occ; d_occ]);
    mon = mmon([t_occ; x_occ; w_occ], degree);
    
    if discrete
        pushforward = subs(mon, vars_new, [f_occ; w_occ]);
        Ay = mom(pushforward - mon);
    else
        Ay = mom(diff(mon, x_occ)*f_occ);
        if ~isempty(var.t)
            Ay = Ay + mom(diff(mon, t_occ));
        end
    end
    X_occ = subs(X, var.x, x_occ);

end

