%test out arguments of peak_estimate routine
%clear
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow
SOLVE = 1;
SAMPLE = 0;
PLOT = 0;

if SOLVE

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('t', 1, 1);
mpol('x', 2, 1);
mpol('d', 1, 1);

dmax = 0.15;
draw = dmax * d;
%support
Xsupp = [];

%dynamics
f = [x(2); -x(1)*(1+ draw) + (1/3).* x(1).^3 - x(2)*(1+ draw)];
% f = [x(2); -x(1)*(1+ draw) + (1/3).* x(1).^3 - x(2)];
X = [];


%initial set
%C0 = [1.2; 0];
C0 = [1.5; 0];
R0 = 0.4;

%C0 = [0.8; 0];
%R0 = 0.1;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%unsafe set

%circle
Cu = [0; -0.5];
%Cu = [2.5; 0];
Ru = 0.5;
c1f = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;


%line
%tilt angle
%order 5: 

theta_c = 5*pi/4;       %p* = -0.1417, beta = [0, 1]
% theta_c = 3*pi/4;       %p* = 0.1935, beta = [0.712, 0.287]
% theta_c = 0;
% theta_c = 3*pi/2;    
% theta_c = 7*pi/4;

w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

%objective to maximize
objective = [c1f; c2f];
%objective = -x(2) - x(1);
%
p_opt = peak_options;
p_opt.var.t = t;
p_opt.var.x = x;
p_opt.var.d = d;

p_opt.disturb = (d^2 <= 1);
p_opt.state_supp = Xsupp;
p_opt.state_init = X0;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;
p_opt.dynamics.discrete = 0;

Tmax_sim = 5;
p_opt.Tmax = Tmax_sim;

p_opt.box = 3;
p_opt.scale = 0;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 4;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val
end

%% now do plots
%there is something wrong with the definition of v. why?
if SAMPLE
rng(50, 'twister')
x0 = C0;

% Nsample = 300;
% Nsample = 150;
Nsample = 20;

s_opt = sampler_options;
s_opt.sample.x = @() sphere_sample(1, 2)'*R0 + C0;
s_opt.sample.d = @() dmax * (2*rand() -1 );
s_opt.Tmax = Tmax_sim;
s_opt.Nd = 1;
s_opt.parallel = 1;

s_opt.mu = 0.2;
tic
out_sim = sampler(out.dynamics, Nsample, s_opt);
sample_time = toc;




end

if PLOT
    %quiver(x, y, xdot, ydot)
    
    nplot = nonneg_plot(out, out_sim);

    %initial and unsafe sets
    theta = linspace(0, 2*pi, 100);
    circ = [cos(theta); sin(theta)];
    %half_theta = linspace(pi/2, 3*pi/2, 100);
    %half_circ = [cos(half_theta); sin(half_theta)];
    
    %initial set    
    X0 = C0 + circ*R0;
    
    %sample from X0 
    rng(33, 'twister')

    Xsample = C0 + circle_sample(Nsample)' * R0;
            
    if theta_c > pi
        theta_c = theta_c - 2 * pi;
    end
    
    
    theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
    circ_half = [cos(theta_half_range); sin(theta_half_range)];
    Xu = Cu + circ_half* Ru;
    
    %expanded set Xu   
    radius_out = sqrt(-peak_val + Ru^2);
    
    phi = asin(-peak_val/radius_out);
    
    th_new_top = phi + theta_c + pi/2; 
    th_new_bot = pi - phi + theta_c + pi/2;
    
    
    if th_new_top > pi
        th_new_top = th_new_top - 2*pi;
        th_new_bot = th_new_bot - 2*pi;
    end
    
    theta_new_range = linspace(th_new_bot, th_new_top + 2*pi, 200);
    circ_margin = [cos(theta_new_range); sin(theta_new_range)];
    Xmargin = Cu + circ_margin*radius_out;
    
    
    
    Xmargin_corner = Cu + circ_margin([1,2],[1,200])*radius_out;
    
    figure(1)
    clf
    hold on
    
    
    for i = 1:Nsample
        if i == 1
            plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
        else
            plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
        end
    end
    
    plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3)
    patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
    
    
    plot(Xmargin(1, :), Xmargin(2, :), 'r--', 'Linewidth', 3)
    plot(Xmargin_corner(1, :), Xmargin_corner(2, :), 'r--', 'Linewidth', 3,   'HandleVisibility','off')

    syms y [2 1]
    vy = out.func.vval(0, y, []) - out.peak_val;

%     fimplicit(vy, [xlim, ylim], ':k', 'DisplayName', 'Invariant Set', 'LineWidth', 3);
        
%     legend({'Trajectories', 'Initial Set', 'Unsafe Set', 'Safety Margin', 'Invariant Set'}, 'location', 'northwest')
           
    xlim([-0.75, 2.25])
    ylim([-1.5, 1.5])
    hold off
    axis square
    xlabel('x_1')
    ylabel('x_2')
    
end


%visualize
