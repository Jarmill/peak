function [event_eval, terminal, direction] = box_event(t, x, box)
%BOX_EVENT Event function for switch_sim used when simulating systems
%   return the event function eventFcn. Yes this is tricky
%is always true.

    %stop integrating when the system falls outside box
    %
    event_eval = all(abs(x) <= box);
    terminal = 1;
    direction = 0;
end
