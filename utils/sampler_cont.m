function [out_sim] = sampler_cont(dynamics, x0, w0, opts)
%SAMPLER_CONT Sample a single trajectory in continuous time. Includes 
%switching, time-independent (parameter) uncertainty, and time-dependent 
%(box and general) uncertainties. Trajectories are simulated until time
%Tmax in the options array. Time-dependent uncertainties change behaviour
%according to an exponential distribution with mean mu.
%
% Input:
%   dynamics:   A struct with fields:
%       f:          System dynamics as a function f: (t, x) -> x' = f(t,x)
%       X:          Support set in space (default everywhere)
%       Tmin:       Minimum time where system is active (default 0)
%       Tmax:       Maximum time where system is active (default Inf)
%   x0:         Initial condition in state
%   w0:         Initial condition for time-independent parametric
%               uncertainty (default is [])
%   opts:       Options for sampler
%
% Output:
%   out_sim:    Data structure with fields
%       t:          time
%       x:          state
%       w:          parametric uncertainty
%       d:          general time-dependent uncertainty
%       b:          box [0,1] time-depedent uncertainty
%       break_sys:  active system in time break
%       break_time: time breaks for system

%gather information about the system

%time range of states
if iscell(dynamics.f)
    nsys = length(dynamics.f);    
else
    nsys = 1;
    dynamics.f = {dynamics.f};    
    dynamics.event = @(t,x) deal(1, 1, 0);
end

%main solving loop
x0_curr = x0;
time_accum = [];
x_accum = [];
time_index= [];

d_accum = [];
b_accum = [];

time_total = 0;

system_choice = [];
time_breaks = 0;

%options = odeset('Events',@events,'OutputFcn',@odeplot,'OutputSel',1,...
%   'Refine',refine);
k = 1;
while time_total < opts.Tmax
    
    %choose a possible system that is admissible for current time/state
    %possible_sys = find(cellfun(@(e) e(time_total, x0_curr), dynamics.event));
    possible_sys = [];
    for i = 1:nsys        
        if opts.Nw > 0
            event_value = dynamics.event{i}(time_total, x0_curr, w0);
        else
            event_value = dynamics.event{i}(time_total, x0_curr);
        end
        if event_value == 1
            possible_sys = [possible_sys; i];
        end
    end
    
    N_possible = length(possible_sys);
    if N_possible == 0
        break
    end
    
    curr_sys_ind = randi([1, N_possible], 1, 1);
    
    
    %This system should be tracked for time time_track or if event is false
    if (nsys == 1) && (opts.Nb == 0) && (opts.Nd == 0)
        time_track = opts.Tmax;
    else
        time_track = exprnd(opts.mu, 1, 1);
    end
    
    curr_sys = possible_sys(curr_sys_ind);
    d_curr = opts.sample.d();
    b_curr = rand(opts.Nb, 1);
    
    curr_f = @(t, x) dynamics.f{curr_sys}(t, x, w0, d_curr, b_curr);
    %currently w is not supported in event function. change this
    curr_event = @(t, x, w) dynamics.event{curr_sys}(t + time_total, x);


    
    time_track_trunc = min(time_track, opts.Tmax - time_total);
    
    %simulate the current system
    curr_ode_options = odeset('Events',curr_event);
    
    [time_curr, x_curr] = opts.odefcn(curr_f, [0, time_track_trunc], x0_curr, curr_ode_options);
    
    
    
    %[time_curr,  x_curr] = ode45(sys_curr, [0, time_range(i)], x0_curr);
    
    %save current trajectory
    %check indices/dimensions
    x0_curr = x_curr(end, :)';
    x_accum = [x_accum; x_curr];
    time_accum = [time_accum; time_curr + time_total];
        
    time_total = time_total + time_curr(end);
    
    %time dependent uncertainty
    if isempty(d_curr)
        d_accum = [d_accum; zeros(size(x_curr, 1), 0)];
    else
%         d_accum = [d_accum; ones(size(x_curr, 1), 1)* d_curr'];  %general
        d_accum = [d_accum; d_curr'];  %general
    end
    b_accum = [b_accum; ones(size(x_curr, 1), 1)* b_curr'];  %box
    system_choice = [system_choice; curr_sys];  %system switching
    time_breaks = [time_breaks; time_total];
    
    time_index= [time_index; k*ones(length(time_curr), 1)];
    k = k + 1;
end

%package up the output
%out = struct('t', times, 'x', x, 'break_sys', sys_id, 'break_time', time_break);
out_sim = struct;

%trajectories
out_sim.t = time_accum;
out_sim.x = x_accum;
out_sim.w = w0;
out_sim.d = d_accum;
out_sim.b = b_accum;

%switching and statistics
out_sim.break_sys = system_choice;
out_sim.break_time = time_breaks;
out_sim.break_index = time_index;
out_sim.Tmax = opts.Tmax;
out_sim.cost = dynamics.cost(x_accum');
%Evaluate (hopefully) nonnegative functions along trajectories
if isfield(dynamics, 'nonneg')    
    out_sim.nonneg = dynamics.nonneg(time_accum', x_accum', w0, d_accum');
end

end

