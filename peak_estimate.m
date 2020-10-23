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

%% Initial and Peak Variables
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
    W0 = subs(Wp, wp, w0);
    W = Wp;
    %Wp = W0 since w is time-independent
    %enforce it explicitly for the sake of compactness
    W = [Wp; W0];  
    
else
    W = [];    
    Wp = [];
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

%% Initial and Peak Measures

%peak time
%mpol('tp', 1);
mup = meas([tp; xp; wp]);
monp = mmon([tp; xp; wp], degree);
yp = mom(monp);

%initial time             
mu0 = meas([t0; x0; w0]); 
mon0 = mmon([t0; x0; w0], degree);
y0 = mom(mon0);


mu_occ = cell(nsys, 1);
mu_occ_sum = 0;
mon_occ = cell(nsys, 1);
v = cell(nsys, 1);
Ay = 0;
%X = cell(nsys, 1);
X_occ = []; %support

%% occupation measure variables

%state
mpol('x_occ', nx, nsys);

%time
if TIME_INDEP
    t_occ = zeros(0, nsys);
else
    mpol('t_occ', 1, nsys);
end


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

%% Box variables

BOX_VAR = nb > 0 && ~options.dynamics.discrete;

if BOX_VAR       
    %state variable x
    xb = cell(nb, 1); %box variables x
    xc = cell(nb, 1); %complement box variables x
    for k = 1:nb
        %define variables
        mpol(strcat('xb', num2str(k)), nx, nsys);
        mpol(strcat('xc', num2str(k)), nx, nsys);
        
        %add variables to cells (janky code but it works)
        xb{k} = eval(strcat('xb', num2str(k)));
        xc{k} = eval(strcat('xc', num2str(k)));
    end
    
    %time variable t
    tb = cell(nb, 1); %box variables x
    tc = cell(nb, 1); %complement box variables x
    if ~TIME_INDEP
        for k = 1:nb
            %define variables
            mpol(strcat('tb', num2str(k)), 1, nsys);
            mpol(strcat('tc', num2str(k)), 1, nsys);

            %add variables to cells (janky code but it works)
            tb{k} = eval(strcat('tb', num2str(k)));
            tc{k} = eval(strcat('tc', num2str(k)));
        end
    end
    
    
    %parameter variable w (time independent)
    wb = cell(nw, 1); %box variables x
    wc = cell(nw, 1); %complement box variables x
    if nw > 0        
        for k = 1:nb
            %define variables
            mpol(strcat('wb', num2str(k)), nw, nsys);
            mpol(strcat('wc', num2str(k)), nw, nsys);

            %add variables to cells (janky code but it works)
            wb{k} = eval(strcat('wb', num2str(k)));
            wc{k} = eval(strcat('wc', num2str(k)));
        end
    end
    
    
    %disturbance variable d (time varying)
    db = cell(dw, 1); %box variables x
    dc = cell(dw, 1); %complement box variables x
    if nd > 0        
        for k = 1:nb
            %define variables
            mpol(strcat('db', num2str(k)), nw, nsys);
            mpol(strcat('dc', num2str(k)), nw, nsys);

            %add variables to cells (janky code but it works)
            db{k} = eval(strcat('db', num2str(k)));
            dc{k} = eval(strcat('dc', num2str(k)));
        end
    end
    
    %package up the variables
    box_var = cell(nb, 1);
    for k = 1:nb
        box_var = struct('tb', tb{k}, 'tc', tc{k}, 'xb', xb{k}, 'xc', xc{k}, ...
            'wb', wb{k}, 'wc', wc{k}, 'db', db{k}, 'dc', dc{k});
    end
end


%this will be tricky. deal with it later.
%involves box control-occupation measures and absolute continuity
%constraints.

%% Form the occupation measures
%support of occupation measure
supp_occ = struct('T', [], 'X', [], 'W', [], 'D', []);

%absolute continuity constraint for box program
abscont = [];

 for i = 1:nsys
    
    %support the valid time range
    if TIME_INDEP
        t_curr = [];
        t_cons = [];
    else            
        t_curr = t_occ(i);
        Tmin_curr = options.dynamics.Tmin(i);
        Tmax_curr = options.dynamics.Tmax(i);
        t_cons = (tp - Tmin_curr/options.Tmax)*(Tmax_curr/options.Tmax - tp) >= 0;            
    end
    
    var_new = struct('t', t_curr, 'x', x_occ(:, i), 'w', w_occ(:, i),...
        'd', d_occ(:, i));

    supp_curr = struct('T', t_cons, 'X', X{i}, 'W', Wp, 'D', options.disturb);
    
    %TODO: finish this part
    if BOX_VAR
        %control-occupation measures with box variables b in [0, 1]
        %efficient formulation only valid if system is continuous
    else        
        [mu_occ{i}, mon_occ{i}, supp_occ_curr, Ay_curr] = occupation_measure(f{i}, ...
            supp_curr, options.var, var_new, degree, options.dynamics.discrete);
    end

    %append supports of current system
    supp_occ.T = [supp_occ.T; supp_occ_curr.T];
    supp_occ.X = [supp_occ.X; supp_occ_curr.X];
    supp_occ.W = [supp_occ.W; supp_occ_curr.W];
    supp_occ.D = [supp_occ.D; supp_occ_curr.D];
    
%     X_occ = [X_occ; t_cons; X_occ_curr];
    Ay = Ay + Ay_curr;
    
    %something about absolute continuity
end

%% Form Constraints and Solve Problem

supp_con = [Xp; X0; supp_occ.X; supp_occ.T; ...
    W; supp_occ.W; supp_occ.D];

%careful of monic substitutions ruining dual variables
Liou = Ay + (y0 - yp);
mom_con = [mass(mu0) == 1; Liou == 0; abscont];

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

out.var = struct('t', tp, 'x', xp, 'w', wp, 'd', options.var.d);


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

%evaluation of beta' p(x)
out.func.cost_beta = @(x) (beta'*cell2mat(cellfun(@(c) c(x), ...
    out.func.cost_all, 'UniformOutput', false)));

if nobj > 1
    out.func.cost = @(x) min(eval(options.obj, xp, x));
else
    out.func.cost = @(x) (eval(options.obj, xp, x));
end

%TODO: missing the 'beta' weights for cost function in this representation
%under multiple cost functions. Check that out later, proper dual
%representation and verification of nonnegativity


%b isn't needed with proper work on control-occupation

if TIME_INDEP
    fval_curr_all = @(t, x, w, d, b) eval(options.dynamics.f{i}, ...
        [xp; wp; options.var.d; options.var.b], [x; w; d; b]);
else
    fval_curr_all = @(t,x, w, d, b) eval(options.dynamics.f{i}, ...
        [tp; xp; wp; options.var.d; options.var.b], [t; x; w; d; b]);
end

%simplify this abomination with [t;x;w;d;b]?
if TIME_INDEP
   out.func.vval = @(t, x, w) eval(v, [xp; wp], [x; repmat(w,size(x, 2))]);    %dual v(t,x,w)
   out.func.Lvval = @(t, x, w, d) eval(Lv, ...
       [xp; wp; options.var.d],  [x; repmat(w,size(x, 2)); d]);
else    
   out.func.vval = @(t, x, w) eval(v, [tp; xp; wp], [t; x; repmat(w,size(x, 2))]);    %dual v(t,x,w)
   out.func.Lvval = @(t, x, w, d) eval(Lv, ...
       [tp; xp; wp; options.var.d],  [t; x; repmat(w,size(x, 2)); d]);
end

%TODO: add functionality for b (sigma)
out.func.nonneg = @(t, x, w, d) ...
                [out.func.vval(t, x, w) + obj_rec; ...
                out.func.Lvval(t, x, w, d); ...
                -out.func.vval(t, x, w) - out.func.cost_beta(x)];
            

out.dynamics.cost = out.func.cost;
out.dynamics.nonneg = out.func.nonneg;
%% Done!

end

function [mu, mon, supp_occ, Ay] = occupation_measure(f, supp, var, var_new, degree, discrete)
%form the occupation measure
%Input:
%   f:        Dynamics (function)
%   supp:     Support set upon which dynamics take place (all variables)
%   var:      variables of problem [t,x,w,d,b]
%   var_new:  New variables for occupation measure
%   degree:   Degree of measure
%   discrete: True (discrete system) or False (continuous system)
%
%Output:
%   mu:       Occupation Measure
%   mon:      monomoials in Liouville constraint
%   supp_occ: Support of measure in all variables
%   Ay:       Adjoint of lie derivative, Liouville

    %var_all = [var.t; x_occ; var.w];

    if nargin < 6
        discrete = 0;
    end
    
    %new support set in the new variables
    supp_occ = struct('T', [], 'X', [], 'W', [], 'D', []);
    
    x_occ = var_new.x;    
    supp_occ.X = subs(supp.X, var.x, x_occ);
        
    if ~isempty(var_new.t)
        t_occ = var_new.t;
        supp_occ.T = subs(supp.T, var.t, t_occ);
    else
        t_occ = [];
    end
    
    if ~isempty(var_new.w)
        w_occ = var_new.w;
        supp_occ.W = subs(supp.W, var.w, w_occ);
    else
        w_occ = [];
    end
    
    if ~isempty(var_new.d)
        d_occ = var_new.d;
        supp_occ.D = subs(supp.D, var.d, d_occ);
    else
        d_occ = [];
    end
    
    %box variables (var.b) are not used here.
    
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

end

function [mu, sigma, sigma_hat, mon, X_occ, Ay, abscont] = occupation_measure_box(f, X, var, var_new, var_box, degree)
%form the occupation measure when there is control-affine box structure:
%f = f0 + sumk bk fk where each bk is a time varying disturbance in [0, 1]
%
%assumptions: length(var.b) > 0 (nontrivial box dynamics)
%             continuous system (discrete does not decompose, pushforward)
%Input:
%   f:        Dynamics (function)
%   X:        Support set upon which dynamics take place
%   var:      variables of problem [t,x,w,d,b]
%   var_new:  New variables for occupation measure
%   var_box:  New variables for box-occupation measure
%   degree:   Degree of measure
%
%Output:
%   mu:         occupation measure (standard)
%   sigma:      box-occupation measure
%   sigma_hat:  complement of box-occupation measure
%   mon:        monomoials in Liouville constraint
%   X:          Support of measure in state
%   Ay:         Adjoint of lie derivative, Liouville (moment constraint)
%   abscont:    Absolute continuity constraint (moment constraint)

    %new measures
    Nb = length(var.b);
    sigma = cell(Nb, 1);
    sigma_c = cell(Nb, 1);
        
    f0 = subs(f, var.b, zeros(length(var.b), 1));
    
    [mu, mon, X_occ, Ay] = occupation_measure(f0, X, var, var_new, degree, 0);

    
   
end