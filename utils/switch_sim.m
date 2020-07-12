function [out] = switch_sim(systems, x0, Tmax, mu)
%SWITCH_SIM Simulate a system that switches between different dynamics from
%time 0 to Tmax. Switches are modeled by an exponential distribution with
%mean mu. This assumes that all dynamics have the same valid region, will
%need to work on the closed cover further later.
%
% Input:
%   systems:    A cell full of function_handles that are the dynamics
%   x0:         Initial condition
%   Tmax:       Maximum time range of simulation
%   mu:         Mean time for system switching
%
% Output:
%   out:        Data structure with fields
%       t:          time
%       x:          state
%       break_sys:  active system in time break
%       break_time: time breaks for system

if nargin < 4
    mu = 1;
end

%gather information about the system
Nsys = length(systems);
n = length(x0);

[sys_id, time_break] = switch_breaks(Nsys, Tmax, mu);
time_range = diff(time_break);
Nbreak = length(sys_id);

%main solving loop
x0_curr = x0;
time_accum = [];
x_accum = [];
for i = 1:Nbreak
    %simulate current system
    sys_curr = systems{sys_id(i)};    
    [time_curr,  x_curr] = ode45(sys_curr, [0, time_range(i)], x0_curr);
    
    %save current trajectory
    %check indices/dimensions
    x0_curr = x_curr(end, :);
    x_accum = [x_accum; x_curr];
    time_accum = [time_accum; time_curr + time_break(i)];
end

%package up the output
%out = struct('t', times, 'x', x, 'break_sys', sys_id, 'break_time', time_break);
out = struct;
out.t = time_accum;
out.x = x_accum;
out.break_sys = sys_id;
out.break_time = time_break;

end

