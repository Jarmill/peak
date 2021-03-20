%find the minimum distance between trajectories of the 'flow' system and a
%half-circle unsafe set

rng(343, 'twister');



SOLVE = 1;
SAMPLE = 0;
PLOT = 1;

n = 2;
order = 4;
d = 2*order;


%% problem parameters
f_func = @(x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
Tmax = 5;
BOX = 3;

if SOLVE
mset clear
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));

mpol('x', 2, 1);
% f = Tmax * [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
f = Tmax*f_func(x);
%initial set
%C0 = [1.2; 0];
C0 = [1.5; 0];
R0 = 0.4;
X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%half-circle set

% Cu = [0; -0.5];
Cu = [0; -0.75];
%Cu = [2.5; 0];
Ru = 0.5;
c1f = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;

theta_c = 5*pi/4;       %p* = -0.1417, beta = [0, 1]
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 




%% set up measures
mpol('t0', 1, 1);
mpol('x0', 2, 1);
mu0 = meas([t0; x0]);

mpol('t_occ', 1, 1);
mpol('x_occ', 2, 1);
mu_occ = meas([t_occ; x_occ]);

mpol('tp', 1, 1);
mpol('xp', 2, 1);
mup = meas([tp; xp]);
%wasserstein
mpol('xw', 2, 1)
mpol('xu', 2, 1)
eta = meas([xw; xu]);

u_cons = subs_vars([c1f; c2f], x, xu);

X0_con = subs_vars((x(1)-C0(1))^2 + (x(2)-C0(2))^2, x, x0) <= R0^2;
% X0_con = (x0 == C0);

%% support constraints
supp_con = [t0 == 0; t_occ*(1-t_occ)>=0; tp*(1-tp) >=0;
    x_occ.^2 <= BOX^2; xp.^2 <= BOX^2; 
    X0_con;
    xw.^2 <= BOX^2; 
    u_cons >= 0;
    ];



%% moment constraints

%liouville constraint
y0 = mom(mmon([t0; x0], d));
yp = mom(mmon([tp; xp], d));


v  = mmon([t_occ; x_occ], d);
f_occ = subs_vars(f, x, x_occ);
Ay = mom(diff(v, t_occ) + diff(v, x_occ)*f_occ); 
Liou = Ay + (y0 - yp);

%marginal between peak and wass
ypx = mom(mmon(xp, d));
ywx = mom(mmon(xw, d));
Wmarg = ywx-ypx;


mom_con = [Liou == 0; Wmarg == 0; mass(mu0)==1];


%objective
mom_dist = mom(sum((xw-xu).^2));
objective = min(mom_dist);

%% solve problem
%Input LMI moment problem
P = msdp(objective, ...
    mom_con, supp_con);

%solve LMIP moment problem
[status, obj] = msol(P);    

%% analyze solutions

dist_rec = sqrt(double(mom_dist));
disp(['distance bound: ', num2str(dist_rec)])
M0 = double(mmat(mu0));
Mocc = double(mmat(mu_occ));
Mp = double(mmat(mup));
Mw = double(mmat(eta));

M0_1 = M0(1:(n+2), 1:(n+2));
Mp_1 = Mp(1:(n+2), 1:(n+2));
Mw_1 = Mw(1:(2*n+1), 1:(2*n+1));

rankp = rank(Mp_1, 1e-3);
rank0 = rank(M0_1, 1e-3);
rankw = rank(Mw_1, 1e-3);

xu_rec = double(mom(xu));
xp_rec = double(mom(xp));
x0_rec = double(mom(x0));


optimal_pt = all([rankp; rank0; rankw]==1);

end 

%% Sample trajectories
if SAMPLE
    Nsample = 20;
    Tmax_sim = 3;
%     sampler = @() circle_sample(1)'*R0 + C0;

    flow_event = @(t, x) box_event(t, x, BOX);
    sample_x = @() circle_sample(1)'*R0 + C0;    

    out_sim = cell(Nsample, 1);
    
    for i = 1:Nsample
        x0_curr = sample_x();
        curr_ode_options = odeset('Events',flow_event);
        [time_curr, x_curr] = ode15s(@(t, x) f_func(x), [0, Tmax], x0_curr, curr_ode_options);
        out_sim{i} = struct('t', time_curr, 'x', x_curr);
        
        %figure out coordinate transformation to evaluate true distances
    end
%     s_opt = sampler_options;
%     s_opt.sample.x = @() circle_sample(1)'*R0 + C0;    
%     s_opt.Tmax = Tmax;

%     s_opt.parallel = 0;
%     s_opt.mu = 0.4;

%     out_sim = sampler(dynamics, Nsample, s_opt);
%     out.optimal = 0;

end

%% Plot flow output
if PLOT
    
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
    
    
    
xu_rec = double(mom(xu));
xp_rec = double(mom(xp));
x0_rec = double(mom(x0));
    scatter(x0_rec(1), x0_rec(2), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);        
    scatter(xp_rec(1), xp_rec(2), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);        
    scatter(xu_rec(1), xu_rec(2), 200, 'sb', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);        
         

end