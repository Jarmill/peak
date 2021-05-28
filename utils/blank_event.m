function [event_eval, terminal, direction] = blank_event(t, x, w)
%BLANK_EVENT

    %stop integrating when the system falls outside support
    %
    event_eval = 1;
    terminal = 1;
    direction = 0;
end