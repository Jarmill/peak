function [out] = peak_estimate(options, order)
%PEAK_ESTIMATE  Find the maximum value of a function p(x) along all
%trajectories starting from an initial set X0 along polynomial dynamics
%f(t, x). This is for finite time [0, Tmax], but if f(t, x) is autonomous
%then the problem can be solved for infinite time.
%   Accomodates switching as well.
%
%Input: options structure (peak_options) with fields:
%   var:        Structure of symbolic variables (@mpol)
%               only the state is nonempty by default
%       t:      time
%       x:      state
%       w:      time-independent parametric uncertainty
%       d:      time-dependent general uncertainty
%       b:      time-dependent box uncertainty (in [0, 1])
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
nvar = 0;

%% Process option data structure
mset clearmeas

%continuous or discrete
if ~isfield(options.dynamics, 'discrete')
    options.dynamics.discrete = 0;
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

%state variable
nx = length(options.var.x);
nvar = nvar + nx;

%parameter variables
nw = length(options.var.w);
nvar = nvar + nw;


nd = length(options.var.d);
nvar = nvar + nd;

nb = length(options.var.b);
% nvar = nvar + nb;

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
    Xp = [subs(options.state_supp, xp, xp_scale); XR]; %that might be the bug (hopefully)
else
    
    %this is where dual_rec breaks
%     Xp = [Xsupp; XR];
    %adding equality constraints here breaks the invariant function
    %dual_rec'monomials. Since X0 and X_occ are constrained within Xsupp,
    %this should not be a problem. If dynamics stay within Xsupp, then Xp
    %should be supported in Xsupp at optimality.
    Xp = [Xsupp; XR];
end
%deal with hanging variables and measures by letting the original x be the
%peak measure

%replace with time breaks
if TIME_INDEP
    t0 = [];
    T0 = [];
else
    mpol('t0', 1, 1);
    T0 = [t0 == 0];
end

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

%time
if TIME_INDEP
    t_occ = zeros(0, nsys);
else
    mpol('t_occ', 1, nsys);
end

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

%% Box variables

%if I was halfway intelligent, I would wrap the occupation measures into a
%class, and use class methods

BOX_VAR = nb > 0 && ~options.dynamics.discrete;

%switch the for loops to index by subsystem 
%for i = 1:nsys
%   xb{i}  = mpol('xb',nx, nb)

if BOX_VAR       
    
    %time variable t
    tb = cell(nsys, 1); %box variables x
    tc = cell(nsys, 1); %complement box variables x
    if ~TIME_INDEP
        for i = 1:nsys
            %define variables
            mpol(strcat('tb', num2str(i)), 1, nb);
            mpol(strcat('tc', num2str(i)), 1, nb);

            %add variables to cells (janky code but it works)
            tb{i} = eval(strcat('tb', num2str(i)));
            tc{i} = eval(strcat('tc', num2str(i)));
        end
    end
    
    %state variable x
    xb = cell(nsys, 1); %box variables x
    xc = cell(nsys, 1); %complement box variables x
    for i = 1:nsys
        %define variables
        mpol(strcat('xb', num2str(i)), nx, nb);
        mpol(strcat('xc', num2str(i)), nx, nb);
        
        %add variables to cells (janky code but it works)
        xb{i} = eval(strcat('xb', num2str(i)));
        xc{i} = eval(strcat('xc', num2str(i)));
    end   
    
    
    %parameter variable w (time independent)
    wb = cell(nsys, 1); %box variables x
    wc = cell(nsys, 1); %complement box variables x
    if nw > 0        
        for i = 1:nsys
            %define variables
            mpol(strcat('wb', num2str(i)), nw, nb);
            mpol(strcat('wc', num2str(i)), nw, nb);

            %add variables to cells (janky code but it works)
            wb{i} = eval(strcat('wb', num2str(i)));
            wc{i} = eval(strcat('wc', num2str(i)));
        end
    end
    
    
    %disturbance variable d (time varying)
    db = cell(nsys, 1); %box variables x
    dc = cell(nsys, 1); %complement box variables x
    if nd > 0        
        for i = 1:nsys
            %define variables
            mpol(strcat('db', num2str(i)), nd, nb);
            mpol(strcat('dc', num2str(i)), nd, nb);

            %add variables to cells (janky code but it works)
            db{i} = eval(strcat('db', num2str(i)));
            dc{i} = eval(strcat('dc', num2str(i)));
        end
    end
    
    %package up the variables
    box_var = cell(nsys, 1);
    for i = 1:nsys
        box_var{i} = struct('tb', tb{i}, 'tc', tc{i}, 'xb', xb{i}, 'xc', xc{i}, ...
            'wb', wb{i}, 'wc', wc{i}, 'db', db{i}, 'dc', dc{i});
    end
    
    sigma_b = cell(nsys, 1);    
    sigma_c = cell(nsys, 1);    
end


%this will be tricky. deal with it later.
%involves box control-occupation measures and absolute continuity
%constraints.

%% Form the occupation measures
%support of occupation measure
supp_occ = struct('T', [], 'X', [], 'W', [], 'D', []);
supp_box = struct('T', [], 'X', [], 'W', [], 'D', []);

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
        [mu_occ{i}, mon_occ{i}, supp_occ_curr, Ay_curr, sigma_b{i}, sigma_c{i}, supp_box_curr, abscont_curr] ...
            = occupation_measure_box(f{i}, supp_curr, options.var, var_new, box_var{i}, degree);
        
        %append absolute continuity constraint
        abscont = [abscont; abscont_curr];
        
        %append supports of current box system
        supp_box.T = [supp_box.T; supp_box_curr.T];
        supp_box.X = [supp_box.X; supp_box_curr.X];
        supp_box.W = [supp_box.W; supp_box_curr.W];
        supp_box.D = [supp_box.D; supp_box_curr.D];
        
    else        
        [mu_occ{i}, mon_occ{i}, supp_occ_curr, Ay_curr] = occupation_measure(f{i}, ...
            supp_curr, options.var, var_new, degree, options.dynamics.discrete);
    end

    %append supports of current system
    supp_occ.T = [supp_occ.T; supp_occ_curr.T];
    supp_occ.X = [supp_occ.X; supp_occ_curr.X];
    supp_occ.W = [supp_occ.W; supp_occ_curr.W];
    supp_occ.D = [supp_occ.D; supp_occ_curr.D];
    
    Ay = Ay + Ay_curr;
    
end

%% Form Constraints and Solve Problem

supp_con = [T0; Tsupp; Xp; X0; supp_occ.T; supp_occ.X;  ...
    W; supp_occ.W; supp_occ.D; ...
    supp_box.T; supp_box.X; supp_box.W; supp_box.D];

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
    
    mom_con = [mom_con; cost_mom_con];
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
[status,obj_rec, m,dual_rec]= msol(P);


%extract moment matrices
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
dual_rec_1 = dual_rec{1};
% dual_rec_v = dual_rec_v((end - length(monp_unscale)+1):end);

n_mon = length(monp_unscale);

dual_rec_v = dual_rec_1(1:n_mon);

if nobj == 1
    beta = 1;
else
    %this is an artifact of gloptipoly not having standard scalar
    %variables. A degree-1 measure is required.
    
    %when 'b' (box uncertainty) is present, check to make sure this is
    %indexed correctly
    beta = dual_rec{2};
end

v = dual_rec_v'*monp_unscale;
if BOX_VAR
    %check with equality constraints of support. There may be some extra
    %dual variables.
    
    %multiplier functions zeta indexed by [box var, subsystem]
    dual_rec_zeta = reshape(dual_rec_1(n_mon + (1:n_mon*nb*nsys)), ...
                            [n_mon, nb, nsys]);
    mon_z = mmon([tp; xp; wp; options.var.d], degree);
    zeta = zeros(nb, nsys)*xp(1); %allocate mpol array
    zeta_sum = zeros(nb, nsys)*xp(1);
    
    %compute Lie derivatives
    Lv =  zeros(1, nsys)*xp(1);
    Lgv = zeros(nb, nsys)*xp(1);
    I = eye(nb);
    for i = 1:nsys
        fi = options.dynamics.f{i};
        
        %base lie derivative (no disturbance)
        fi0 = subs(fi, options.var.b, zeros(nb, 1));
        
        Lv(i) = diff(v, xp)*fi0;
        if ~TIME_INDEP
            Lv(i) = Lv(i) + diff(v, tp);
        end
        
        %box lie derivative (for each box input channel)
        for k = 1:nb
            zeta(k, i) = dual_rec_zeta(:, k, i)'*mon_z;
            zeta_sum(i) = zeta_sum(i) + zeta(k, i);
            
            fik = subs(fi, options.var.b, I(k, :)) - fi0;
            
            Lgv(i, k) = diff(v, xp) * fik;                        
        end
    end
    
    %formulate polynomials that should be nonnegative
    mu_con = (Lv - zeta_sum)';
    sigma_con = reshape(Lgv + zeta, [], 1);
    sigmac_con = reshape(zeta, [], 1);
    
    boxcon = [mu_con; sigma_con; sigmac_con];       
else
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


end


%TODO: process Lv into the box version with 'zeta' from dual_rec

%% Output results to data structure
%out = 1;

out = struct;

%recover optima
out.order = order;
out.peak_val = obj_rec;
out.optimal = (rank0 == 1) && (rankp == 1);
out.x0 = x0_rec;
out.xp = xp_rec;

%moment matrices
out.M0 = M0;
out.Mp = Mp;

%occupation measure
out.Mocc=cellfun(@(m) double(mmat(m)), mu_occ, 'UniformOutput', false);

%box occupation measures
if BOX_VAR
    %vertcat(cell{:}) flattens the nested cell array
    out.Mbox =cellfun(@(m) double(mmat(m)), vertcat(sigma_b{:}), 'UniformOutput', false);
    out.Mboxc =cellfun(@(m) double(mmat(m)), vertcat(sigma_c{:}), 'UniformOutput', false);
end


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

%% Dynamics and Event functions
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
    
    
    event_curr = @(t,x,w) support_event(t, x, Xval_curr, ...
            options.dynamics.Tmin(i), options.dynamics.Tmax(i)); %should w be present here?
    
%     out.func.fval_all{i} = fval_curr_all;
    out.func.fval{i} = fval_curr_all;
    out.func.Xval{i} = Xval_curr;        
    out.func.event{i} = event_curr;
end

out.dynamics = struct;
out.dynamics.f = out.func.fval;
% out.dynamics.f_all = out.func.fval_all;
out.dynamics.event = out.func.event;
out.dynamics.discrete = options.dynamics.discrete;
out.dynamics.time_indep = TIME_INDEP;
% event_all = @(tt, xt) cell2mat(cellfun(@(ex) ex(tt, xt), out.func.event, 'UniformOutput', false));

%% Dual variable functions (value, nonnegativity)
out.func = struct;
out.func.dual_rec = dual_rec;
out.func.v = v;
out.func.beta = beta;
if BOX_VAR
    out.func.Lv = Lv;
    out.func.Lgv = Lgv;
    out.func.zeta = zeta;    
else
    out.func.Lv = Lv;
end

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

%simplify this abomination with [t;x;w;d]?
if TIME_INDEP
   out.func.vval = @(t, x, w) eval(v, [xp; wp], [x; repmat(w,size(x, 2))]);    %dual v(t,x,w)
   out.func.Lvval = @(t, x, w, d) eval(Lv, ...
       [xp; wp; options.var.d],  [x; repmat(w,size(x, 2)); d]);
else    
   out.func.vval = @(t, x, w) eval(v, [tp; xp; wp], [t; x; repmat(w,size(x, 2))]);    %dual v(t,x,w)
   out.func.Lvval = @(t, x, w, d) eval(Lv, ...
       [tp; xp; wp; options.var.d],  [t; x; repmat(w,size(x, 2)); d]);
end
if BOX_VAR
    if TIME_INDEP
       out.func.boxcon = @(t, x, w, d) eval(boxcon, ...
           [xp; wp; options.var.d],  [x; repmat(w,size(x, 2)); d]);
    else           
       out.func.boxcon=   @(t, x, w, d) eval(boxcon, ...
           [tp; xp; wp; options.var.d],  [t; x; repmat(w,size(x, 2)); d]);
    end
end

%TODO: add functionality for b (zeta)
if BOX_VAR
    out.func.nonneg = @(t, x, w, d) ...
                    [out.func.vval(t, x, w) + obj_rec; ...                                        
                    out.func.boxcon(t, x, w, d); ...
                    -out.func.vval(t, x, w) - out.func.cost_beta(x)];
else
    out.func.nonneg = @(t, x, w, d) ...
                    [out.func.vval(t, x, w) + obj_rec; ...
                    out.func.Lvval(t, x, w, d); ...
                    -out.func.vval(t, x, w) - out.func.cost_beta(x)];
end        

            
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
    vars_sub = [t_occ; x_occ; w_occ];
    
    f_occ = subs(f, vars_prev, vars_new);
    
    mu = meas([t_occ; x_occ; w_occ; d_occ]);
    mon = mmon([t_occ; x_occ; w_occ], degree);
    
    if discrete
        pushforward = subs(mon, vars_sub, [f_occ; w_occ]);
        Ay = mom(pushforward - mon);
    else
        Ay = mom(diff(mon, x_occ)*f_occ);
        if ~isempty(var.t)
            Ay = Ay + mom(diff(mon, t_occ));
        end
    end    

end


function [mu, mon_occ, supp_occ, Ay, sigma, sigma_c, supp_box, abscont] = ...
    occupation_measure_box(f, supp, var, var_new, var_box, degree)
%form the occupation measure when there is control-affine box structure:
%f = f0 + sumk bk fk where each bk is a time varying disturbance in [0, 1]
%
%assumptions: length(var.b) > 0 (nontrivial box dynamics)
%             continuous system (discrete does not decompose, pushforward)
%Input:
%   f:        Dynamics (function), disturbance-affine in box variable b
%   supp:     Support set upon which dynamics take place
%   var:      variables of problem [t,x,w,d,b]
%   var_new:  New variables for occupation measure [t,x,w,d]
%   var_box:  New variables for box-occupation measure [t,x,w,d] for sigma
%             and sigma_c (standard and complement disturbance-occupation)
%   degree:   Degree of measure
%
%Output:
%   mu:         occupation measure (standard)
%   mon_occ:    monomoials in Liouville constraint
%   supp:       Support of occupation measure variables
%   Ay:         Adjoint of lie derivative, Liouville (moment constraint)
%   sigma:      box-occupation measures and their complements
%   supp_box:   Support of box occupation variables
%   abscont:    Absolute continuity constraint (moment constraint)

    %new support set in the new variables
    supp_box = struct('T', [], 'X', [], 'W', [], 'D', []);

    %new measures
    nb = length(var.b);
    sigma = cell(nb, 1);
    sigma_c = cell(nb, 1);
        
    %base occupation measure (with no box disturbance)
    f0 = subs(f, var.b, zeros(nb, 1));
    
    [mu, mon_occ, supp_occ, Ay] = occupation_measure(f0, supp, var, var_new, degree, 0);
    y_occ = mom(mon_occ); 
    
    %absolute continuity constraint of sigmak with respect to mu
    abscont = [];
    
    %process the box variables
    
    
    %loop through box variables
    I = eye(nb);
    
    for k = 1:nb
        %variables and measures of disturbance occupation
        
        %t time
        if ~isempty(var_new.t)
            tbk = var_box.tb(k);
            tck = var_box.tc(k);
            supp_box.T = [supp_box.T; subs(supp.T, var.t, tbk); subs(supp.T, var.t, tck)];
        else
            tbk = [];
            tck = [];
        end
        
        %x state
        xbk = var_box.xb(:, k);
        xck = var_box.xc(:, k);
        supp_box.X = [supp_box.X; subs(supp.X, var.x, xbk); subs(supp.X, var.x, xck)];

        
        %w parameter
        if ~isempty(var_new.w)
            wbk = var_box.wb(k);
            wck = var_box.wc(k);
            supp_box.W = [supp_box.W; subs(supp.W, var.w, wbk); subs(supp.W, var.w, wck)];
        else
            wbk = [];
            wck = [];
        end

        %d disturbance
        if ~isempty(var_new.d)
            dbk = var_box.db(k);
            dck = var_box.dc(k);
            supp_box.D = [supp_box.D; subs(supp.D, var.d, dbk); subs(supp.D, var.d, dck)];
        else
            dbk = [];
            dck = [];
        end
        
        var_b = [tbk; xbk; wbk; dbk];
        var_c = [tck; xck; wck; dck];
        
        sigma{k}  = meas(var_b);    %disturbance occupation
        sigma_c{k} = meas(var_c);   %complement of disturbance occupation
        
        %monomials and moments
        mon_b = mmon(var_b, degree);
        yb = mom(mon_b);
        yc = mom(mmon(var_c, degree));
        
        %absolute continuity
        %check the sign with regards to the dual variable zeta
        abscont = [abscont; yb + yc - y_occ == 0];
        
        %dynamics on box variable
        fk = subs(f, var.b, I(:, k)) - f0;
        
        %<d/dx_k v(t, x, w), sigma_k >
        Ayk = mom(diff(mon_b, xbk) * fk);
        
        Ay = Ay + Ayk;
    end

end