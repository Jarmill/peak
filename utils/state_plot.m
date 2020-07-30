function [fig] = state_plot(out, out_sim, out_sim_peak)
%STATE_PLOT Plot states of system trajectories, depending on the number of
%states. Passes to specialized functions
%out: information about the recovered solution
%out_sim: function evaluation on random sample trajectories
%out_sim_peak{1}: optional argument, function evaluation if peak is found at
%x0
nx = size(out_sim{1}.x, 2);

if nargin == 2
    out_sim_peak = [];
end

if nx == 2
    fig = state_plot_2(out, out_sim, out_sim_peak);
elseif nx == 3
    fig = state_plot_3(out, out_sim, out_sim_peak);
else
    fig = state_plot_N(out, out_sim, out_sim_peak);
end

end

