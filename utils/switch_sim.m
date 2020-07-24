function [out] = switch_sim(dynamics, x0, Tmax, mu, odefcn)
%SWITCH_SIM Simulate a system that switches between different dynamics from
%time 0 to Tmax. Switches are modeled by an exponential distribution with
%mean mu. This assumes that all dynamics have the same valid region, will
%need to work on the closed cover further later.
%
% Input:
%   dynamics:   A struct with fields:
%       f:          System dynamics as a function f: (t, x) -> x' = f(t,x)
%       X:          Support set in space (default everywhere)
%       Tmin:       Minimum time where system is active (default 0)
%       Tmax:       Maximum time where system is active (default Inf)
%   x0:         Initial condition
%   Tmax:       Maximum time range of simulation
%   mu:         Mean time for system switching
%   odefcn:     Function handle to the ode solver (default ode15s to deal
%               with stiffness)
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

if nargin < 5
    odefcn = @ode15s;
end

%gather information about the system

%time range of states
if iscell(dynamics.f)
    nsys = length(dynamics.f);    
else
    nsys = 1;
    dynamics.f = {dynamics.f};    
    dynamics.event = @(t,x) deal(1, 1, 0);
end

%support in space and time
nx = length(x0);

%nominal time breaks, ignoring dynamics.X{i}
%[sys_id, time_break] = switch_breaks(Nsys, Tmax, mu);
%time_range = diff(time_break);
%Nbreak = length(sys_id);

%main solving loop
x0_curr = x0;
time_accum = [];
x_accum = [];
time_index= [];


time_total = 0;

system_choice = [];
time_breaks = 0;

%options = odeset('Events',@events,'OutputFcn',@odeplot,'OutputSel',1,...
%   'Refine',refine);
k = 1;
while time_total < Tmax
    
    %choose a possible system that is admissible for current time/state
    %possible_sys = find(cellfun(@(e) e(time_total, x0_curr), dynamics.event));
    possible_sys = [];
    for i = 1:nsys
        [event_value, ~, ~] = dynamics.event{i}(time_total, x0_curr);
        if event_value == 1
            possible_sys = [possible_sys; i];
        end
    end
    N_possible = length(possible_sys);
    curr_sys_ind = randi([1, N_possible], 1, 1);
    
    curr_sys = possible_sys(curr_sys_ind);
    curr_f = dynamics.f{curr_sys};
    curr_event = @(t, x) dynamics.event{curr_sys}(t + time_total, x');
    
    %This system should be tracked for time time_track or if event is false
    time_track = exprnd(mu, 1, 1);
    time_track_trunc = min(time_track, Tmax - time_total);
    
    %simulate the current system
    curr_ode_options = odeset('Events',curr_event);
    
    %[time_curr, x_curr] = ode15s(curr_f, [0, time_track_trunc], x0_curr, curr_ode_options);
    %[time_curr, x_curr] = ode45(curr_f, [0, time_track_trunc], x0_curr, curr_ode_options);
    [time_curr, x_curr] = odefcn(curr_f, [0, time_track_trunc], x0_curr, curr_ode_options);
    
    
    
    %[time_curr,  x_curr] = ode45(sys_curr, [0, time_range(i)], x0_curr);
    
    %save current trajectory
    %check indices/dimensions
    x0_curr = x_curr(end, :);
    x_accum = [x_accum; x_curr];
    time_accum = [time_accum; time_curr + time_total];
        
    time_total = time_total + time_curr(end);
    
    system_choice = [system_choice; curr_sys];
    time_breaks = [time_breaks; time_total];
    
    time_index= [time_index; k*ones(length(time_curr), 1)];
    k = k + 1;
end

%package up the output
%out = struct('t', times, 'x', x, 'break_sys', sys_id, 'break_time', time_break);
out = struct;
out.t = time_accum;
out.x = x_accum;
out.break_sys = system_choice;
out.break_time = time_breaks;
out.break_index = time_index;

%TODO: modify for evaluation on time varying system
if isfield(dynamics, 'nonneg')    
    out.nonneg = dynamics.nonneg(x_accum');
end

end

