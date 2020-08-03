%pendulum with unit constants (length, mass, gravity)
rng(300, 'twister')

CONTROLLER = 2;


%lqr realization
A = [0 1; 1 0];
B = [0; -1];

Q = diag([2, 1]);
R = 1;

[K, S, clp] = lqr(A, B, Q, R);


%this is a stable choice when |th-pi|<=pi/4 and |w| <= 1.
%u_lqr = @(x) K*[sin(x(1) - pi); x(2)];
u_lqr = @(x) K*[-sin(x(1)); x(2)];


%swing-up
k = 1; %swing-up gain
u_swing = @(x) -k * x(2) * ((0.5*x(2)^2 - cos(x(1))) - 1);



if CONTROLLER == 0
    %lqr    
    th_max = pi/4;
    w_max = 1;
    u = u_lqr;
elseif CONTROLLER == 1    
    th_max = pi;
    w_max = 1;

    u = u_swing;
else
    th_max = pi;
    w_max = 0;

    eps = 0.1;
    
    
    u = @(x) u_hybrid(x, S, eps, u_lqr, u_swing);
end

%% Simulate
Nsample = 100;
Tmax_sim = 20;

sampler = @() [pi + (2*rand() - 1)* th_max; (2*rand()-1)*w_max];
%out_sim = pend_sampler(sampler, Nsample, Tmax_sim);
out_sim = pend_sampler(sampler, Nsample, Tmax_sim, [], 0, u);

%% Plot
figure(1)
clf
for i = 1:length(out_sim)
    subplot(5, 1, 1)
    if i== 1
        hold on
        title('Cosine of Angle')
        xlabel('time')
        ylabel('cos(\theta)')
    end
    plot(out_sim{i}.t, out_sim{i}.x(:, 1), 'c')

    subplot(5, 1,2)
    if i== 1
        hold on
        title('Sine of Angle')
        xlabel('time')
        ylabel('sin(\theta)')
    end
    plot(out_sim{i}.t, out_sim{i}.x(:, 2), 'c')

    subplot(5, 1,3)
    if i== 1
        hold on
        title('Angular Velocity')
        xlabel('time')
        ylabel('\omega')
    end
    hold on
    plot(out_sim{i}.t, out_sim{i}.x(:, 3), 'c')
    
    out_sim{i}.energy = 0.5*out_sim{i}.x(:, 3).^2 - out_sim{i}.x(:, 1);
    
    subplot(5, 1,4)
    if i== 1
        hold on
        title('Energy')
        xlabel('time')
        ylabel('E')
    end
    hold on
    plot(out_sim{i}.t, out_sim{i}.energy, 'c')
    
    
    out_sim{i}.u = zeros(size(out_sim{i}.t));
    for k = 1:length(out_sim{i}.u)
        out_sim{i}.u(k) = u(trig_process(out_sim{i}.x(k, :)));
    end
    subplot(5, 1,5)
    if i== 1
        hold on
        title('Input')
        xlabel('time')
        ylabel('u')
    end
    hold on
    plot(out_sim{i}.t, out_sim{i}.u, 'c')
end

subplot(5, 1, 1)
plot(xlim, [-1, -1], ':k')
subplot(5, 1, 2)
plot(xlim, [0, 0], ':k')
subplot(5, 1, 3)
plot(xlim, [0, 0], ':k')
subplot(5, 1, 4)
plot(xlim, [1, 1], ':k')
subplot(5, 1, 5)
plot(xlim, [0, 0], ':k')


% % % 
% if (out.optimal == 1)
%     x0_angle = [atan2(out.x0(2), out.x0(1)); out.x0(3)];
%     out_sim_peak = pend_sampler(x0_angle, 1, Tmax_sim, out.dynamics.nonneg, b);
%     nplot = nonneg_plot(out_sim, out_sim_peak);
% % 
%     splot = state_plot_3(out, out_sim, out_sim_peak);
% % %     splot = state_plot_N(out, out_sim, out_sim_peak);
% else
%     nplot = nonneg_plot(out_sim);
% %     splot = state_plot_2(out, out_sim);
%     splot3 = state_plot_3(out, out_sim);
% %     splot = state_plot_N(out, out_sim);
% end
% % 

function x = trig_process(x_trig)
    x = [atan2(x_trig(2), x_trig(1)); x_trig(3)];
end
    
function u_out = u_hybrid(x, S, eps, u_lqr, u_swing)
%hybrid controller.
%if abs(energy -1) >= eps, use u_swing
%if abs(energy-1) <= eps and norm(x)^2 <= delta, use u_lqr
%else, wait
energy_gap = 0.5*x(2)^2 - cos(x(1)) - 1;

if abs(energy_gap) >= eps
    u_out = u_swing(x);
else
    %if norm(x)^2 <= delta
    %sx = [sin(x(1) - pi); x(2)];
    sx = [-sin(x(1)); x(2)];
    Sform = sx'*S*sx;
    if Sform <= 2
        u_out = u_lqr(x);
    else
        u_out = 0;
    end
end


%three states:
%swing-up
%wait
%lqr

end