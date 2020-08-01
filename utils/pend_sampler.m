
function x0 = pend_sampler(th_max, w_max)
    %uniformly sample from an initial pendulum state
    % |th| <= th_max
    % |w|  <= w_max
    th = (rand()*2 - 1)*th_max;
    w = (rand()*2 - 1)*w_max;
    
    x0 = [cos(th); sin(th); w];
end