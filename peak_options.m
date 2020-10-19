classdef peak_options < handle
    %attributes of peak estimation routine
    %a default set of options
    properties
        
        %% properties of run
        %terminal time   
        Tmax(1,1) double{mustBePositive}  = Inf;           
        
        %function dynamics
        %f: dynamics
        %X: space on which dynamics are valid (arbitrary switching)
        dynamics = struct('f', [], 'X', {}, 'Tmin', [], 'Tmax', [], 'discrete', 0)
        
        
        %objective to minimize
        %could be a function (single objective)
        %or a cell of functions (minimum of entries)        
        obj = [];
        
        %iterative cuts
        %output of prior problem may be infeasible (some matrices not PSD)
        %find a feasible point with a cost of approximately prev_cost
        prev_cost = [];
        
        %% Variables and descriptors
        %variables (array of symbolic variables)
%           t:  time
%           x:  state
%           w:  time-independent uncertainty (parameters)
%           d:  time-dependent uncertainty (general)
%           b:  time-dependent uncertainty (constrained to a box)
        var = struct('t', [], 'x', [], 'w', [], 'd', [], 'b', []);
        
        %% support sets
        %type @mpol/supcon
        %(X, T): Xt is the support of trajectories at time Tt
        state_init = []; %initial state at time 0
        state_supp = []; %states to consider
        R = 10; %sum(x.^2) <= R^2
        
        %Coordinate ranges for variables for scaling
        %state variables lie in a box (utils/box_process)
        
        %box is scalar B:           -B      <= x_i <= B
        %box is [Bmin, Bmax]:       -Bmin   <= x_i <= Bmax
        %box is [B_i]:              -Bi     <= x_i <= B_i
        %box is [Bmin_i, Bmax_i]    -Bmin_i <= x_i <= Bmax_i
        box = 1; 
        
        scale = 0; %should variables be scaled to [-1,1] (state) and [0,1] (time)
                        
        param = [];       %parameters  w (@mpol/supcon)  
        disturb = [];     %uncertainty d (@mpol/supcon)  
        
        %not needed, done implicitly.
        %later will add conversion code to handle this
%         disturb_box = []; %uncertainty b ([0,1]^number of b parameters)  
        
        %% additional options
        solver = 'mosek';
        
        %what is the tolerance for identifying a matrix as rank-1?
        rank_tol(1,1) double{mustBePositive} = 1e-3; 
        
        
    end
end