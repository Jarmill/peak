function [system, time_breaks] = switch_breaks(N_system, Tmax, mu)
%SWITCH_BREAKS Switches between N_system different subsystems until time
%Tmax. Switch times governed by exponential distribution with mean mu.
% Input:
%   N_system:   Number of subsystems
%   Tmax:       Time period
%   mu:         Mean time to switch (mu = 1/lambda)
%
% 
% Output:
%   system:     Which system is currently active
%   times:      Time at which system becomes active

%sample time breaks
if N_system == 1
    system = 1;
    time_breaks = [0 Tmax];
else
    step_range = ceil(Tmax/mu);
    times = 0;
    sum_times = 0;

    while sum_times <= Tmax
        times_curr = exprnd(mu, step_range, 1);
        sum_times = sum_times + sum(times_curr);
        times = [times; times_curr];
    end

    %truncate times to within Tmax
    cumu_times = cumsum(times);
    %viol_ind = find(cumu_times > Tmax);
    time_breaks = [cumu_times(cumu_times < Tmax); Tmax];


    N_breaks = length(time_breaks) - 1;

    %now sample subsystems. Initial subsystem and jumps
    sys_init = randi([1, N_system], 1, 1);
    sys_jumps = randi([1, N_system-1], N_breaks-1, 1);

    sys_step = [sys_init; sys_jumps];
    sys_cumu = cumsum(sys_step);

    system = mod(sys_cumu, N_system)+1;
end

end

