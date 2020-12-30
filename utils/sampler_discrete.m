function [out_sim] = sampler_discrete(dynamics, x0, w0, opts)
%SAMPLER_DISCRETE Sample a single trajectory in discrete time. Includes 
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
Tmax = opts.Tmax;
T = (0:Tmax);
X = zeros(Tmax + 1, length(x0));
X(1, :) = x0;
x_curr = x0;
system_choice = [];

% d_accum = [];
d_accum = opts.sample.d();


for k = 2:(Tmax+1)
    
    %find which systems are valid
    possible_sys = [];
    for i = 1:nsys        
        event_value = dynamics.event{i}(k, x_curr, w0);

        if event_value == 1
            possible_sys = [possible_sys; i];
        end
    end
    N_possible = length(possible_sys);
    
    if N_possible == 0
        break
    end
    
    curr_sys_ind = randi([1, N_possible], 1, 1);
    
    %choose which system to use. curr_event not needed, validity evaluated
    %on next state iteration
    curr_sys = possible_sys(curr_sys_ind);            
    system_choice(k-1) = curr_sys;
    
    %sample the time-varying disturbance
    d_curr = opts.sample.d();
    curr_f = @(t, x) dynamics.f{curr_sys}(t, x, w0, d_curr, []);
    
    %take the step and store results
    x_next = curr_f(k, x_curr);
    X(k, :) = x_next;
    if isempty(d_curr)
        d_accum = [d_accum; zeros(size(x_curr, 1), 0)];
    else
        %not sur4e about this line
        d_accum = [d_accum; d_curr'];  %general
    end
    x_curr = x_next;    
end


%pad an additional disturbance
if isempty(d_curr)
    d_accum = [d_accum; zeros(size(x_curr, 1), 0)];
else
    %not sur4e about this line
    d_accum = [d_accum; d_curr'];  %general
end



%package up the output
%out = struct('t', times, 'x', x, 'break_sys', sys_id, 'break_time', time_break);
out_sim = struct;

%trajectories
out_sim.t = T;
out_sim.x = X;
out_sim.w = w0;
out_sim.d = d_accum;
%switches
out_sim.break_sys = system_choice;

out_sim.cost = dynamics.cost(X');

%pad d with zeros in case of early termination (system exits X)
d_accum = [d_accum; zeros(opts.Tmax +1 - size(d_accum, 1), size(d_accum, 2))];

if isfield(dynamics, 'nonneg')    
    out_sim.nonneg = dynamics.nonneg(T', X', w0, d_accum');
end


end

