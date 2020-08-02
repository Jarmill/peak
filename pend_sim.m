function [out_sim] = pend_sim(x0, Tmax, nonneg, b, u, odefcn)
%PEND_SIM Simulate a pendulum with unit constants under control law u
%in terms of theta and theta_dot, not trig representation
%time 0 to Tmax. Switches are modeled by an exponential distribution with
%mean mu. This assumes that all dynamics have the same valid region, will
%need to work on the closed cover further later.
%
% Input:
%   x0:         Initial condition
%   Tmax:       Maximum time range of simulation
%   nonneg:     A function for (hopefully) nonnegative quantities along
%               trajectories
%   b:          friction
%   u:          Control law x -> u(x) 
%   odefcn:     Function handle to the ode solver (default ode15s to deal
%               with stiffness)
%
% Output:
%   out_sim:    Data structure with fields
%       t:          time
%       x:          state

if nargin < 3
    nonneg = [];
end

if nargin < 4
    u = @(x) 0;
end

if nargin < 5
    odefcn = @ode15s;
end

%closed loop dynamics
curr_f = @(t, x) unit_pend_dynamics(t, x, b);
[time_accum, x_accum] = odefcn(curr_f, [0, Tmax], x0);
   

%package up the output
%out = struct('t', times, 'x', x, 'break_sys', sys_id, 'break_time', time_break);
out_sim = struct;
out_sim.t = time_accum;
out_sim.theta = x_accum;
out_sim.Tmax = Tmax;

x_trig = [cos(x_accum(:, 1)), sin(x_accum(:, 1)), x_accum(:, 2)];
out_sim.x = x_trig;

%TODO: modify for evaluation on time varying system
if ~isempty(nonneg)    
    if nargin(nonneg) == 2
        out_sim.nonneg = nonneg(time_accum', x_trig');
    else
        out_sim.nonneg = nonneg(x_trig');
    end
end

end

