%% peak impulse response for a state space system
%from the Chesi papers

%right now deterministic system
%then move to uncertainty

nsys = 3;

SOLVE = 1;
PLOT = 1;
MD = 200; %mesh density for 3d
Tmax = 20;
if SOLVE
if nsys == 1
    %DC motor example
    %constant Jm
    bm = 0.2;
    Kt = 1;
    La = 0.5;
    Ra = 1;
    Ke = 0.5;
    Jm = 1;
    
    A = [0, 1, 0; 0, -bm/Jm, Kt/Jm; 0, -Ke/La, -Ra/La];
    B = [0; 0; 1/La];
    C = [1, 0, 0];
    
    sys = ss(A, B, C, 0);
    
elseif nsys == 2
    %random system
    n = 6;
    sys = rss(n);
else
    %second order system
    %true upper bound is p* = 0.645
    sys = ss([0 1; -0.5 -1], [0; 1], [1 0], 0);
end


%% Gloptipoly setup

%gloptipoly options
order = 2;
n = size(sys.A, 1);
d = 2*order;
rank_tol = 1e-2;


mset clear; warning('off','YALMIP:strict')
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));


%gloptipoly measures
mpol('x', n);  mu  = meas(x);   %occupation measure    
mpol('xp', n); mup = meas(xp); %peak measure

%test functions (monomials)

%v0 = mmon(x0, d);
v  = mmon(x, d);
vp = mmon(xp, d);

%unknown moments of measures
yp = mom(vp);

t = 0;

f = sys.A*x;

%Liouville Equation 
Ay = mom(diff(v, x)*f); 

%initial measure is a point mass at B
powers = genPowGlopti(n, d);
%v0 =  subs(v, x, sys.B);
y0 = prod((sys.B').^powers, 2);

%% Gloptipoly constraints
%measure
%Liouville Equation
Liou = Ay + (y0 - yp);

mom_con = (Liou == 0);
%by Liouville, mass(mup) = mass(mu0) = 1 (atomic)



%support
R = 4*norm(sys.B)^2;
X  = (x'*x <= R^2);
Xp = (xp'*xp <= R^2);

supp_con = [X, Xp];

cost = (sys.C*xp)^2;

objective = max(cost);


%% Solve LMI and recover solution
    
%Input LMI moment problem
P = msdp(objective, ...
    mom_con, supp_con);

%solve LMI moment problem    
[status,obj,m,dual_rec]= msol(P);
[model.A, model.b, model.c, model.K, model.b0, model.s] = msedumi(P);

%extract solutions
%gamma_val = -obj;
peak_val = sqrt(obj);
Mp = double(mmat(mup));
Mp_1 = Mp(1:3, 1:3);
rankp = rank(Mp_1, rank_tol);
xp_out = double(mom(xp));
    
%get the dual polynomial
syms xc [n, 1];
% 
% if n == 2
%     vv = conj(monolistYalToGlop([xc; yc], d));
%     p = dual_rec'*vv;
%     Lp = [diff(p, xc) diff(p, yc)]*sys.A*[xc; yc];
% elseif n==3
%     %n = 3
%     vv = conj(monolistYalToGlop([xc; yc; zc], d));
%     p = dual_rec'*vv;
%     Lp = [diff(p, xc) diff(p, yc) diff(p, zc)]*sys.A*[xc; yc; zc];
% end

vv = conj(monolistYalToGlop(xc, d));
p = dual_rec{1}'*vv;
%Lp = diff(p, xc)*sys.A*xc;
Lp = jacobian(p, xc)*sys.A*xc;

pval = matlabFunction(p);

%no time varying
Lpval = matlabFunction(Lp);

end
%% Plot the results

if PLOT
figure(1)
clf
xtraj = struct;
[xtraj.y, xtraj.t, xtraj.x] = impulse(sys, Tmax);
margin = 0.2;
xmin = min(xtraj.x) - 0.3;
xmax = max(xtraj.x) + 0.3;
trange = xtraj.t([1, end]);

if (n==2) || (n==3)
subplot(1,2,1)
hold on
if n == 3

    plot3(xtraj.x(:, 1), xtraj.x(:, 2), xtraj.x(:, 3), 'DisplayName', 'state');
    scatter3(sys.B(1), sys.B(2), sys.B(3), 100, 'ok', 'DisplayName', 'Initial Condition', 'LineWidth', 2)
    %fi = fimplicit3(p + obj, [xmin(1), xmax(1), xmin(2), xmax(2),xmin(3), xmax(3)], 'MeshDensity',MD, ...
    %        'EdgeColor', 'None','FaceColor', 'k', 'FaceAlpha', 0.3, ...
    %        'DisplayName', 'Safety Contour');
    if rankp == 1
        scatter3(xp_out(1), xp_out(2), xp_out(3), 100, '*k', 'DisplayName', 'Peak Value', 'LineWidth', 2)
    end
    zlabel('x_3')
else
    %n = 2
    plot(xtraj.x(:, 1), xtraj.x(:, 2), 'DisplayName', 'state', 'LineWidth', 3);
    scatter(sys.B(1), sys.B(2), 100, 'ok', 'DisplayName', 'Initial Condition', 'LineWidth', 2)
    
    p_neg = line_range([sys.C, -peak_val], [xmin(1), xmax(1)], [xmin(2), xmax(2)]);
    p_pos = line_range([sys.C, peak_val], [xmin(1), xmax(1)], [xmin(2), xmax(2)]);
    
    if rankp == 1
        scatter(xp_out(1), xp_out(2), 100, '*k', 'DisplayName', 'Peak Value', 'LineWidth', 2)
    end
    if ~isempty(p_neg)
        plot([p_neg{1}(1), p_neg{2}(1)],  [p_neg{1}(2), p_neg{2}(2)], 'r--', 'Linewidth', 3, 'DisplayName', 'Peak Certification')
    end
    if ~isempty(p_pos)
        plot([p_pos{1}(1), p_pos{2}(1)],  [p_pos{1}(2), p_pos{2}(2)], 'r--', 'Linewidth', 3, 'HandleVisibility','off')
    end

    
    fi2 = fimplicit(p  + obj, [xmin(1), xmax(1), xmin(2), xmax(2)], ...
            ':k', 'DisplayName', 'Peak Contour', 'LineWidth', 3, 'MeshDensity', 150);    
end
xlabel('x_1')
ylabel('x_2')
title(['Impulse Response State, order = ', num2str(order)])
axis tight
legend('location', 'northeast')
    

subplot(1,2,2)
end
hold on
title(['Impulse Response Output, peak = ', num2str(peak_val, 3)])
plot(xtraj.t, xtraj.y, 'Linewidth', 3, 'DisplayName', 'Impulse Response')
plot(trange, peak_val*[1,1], 'r--', 'Linewidth', 3, 'DisplayName', 'Peak Certification')
plot(trange, -peak_val*[1,1], 'r--', 'Linewidth', 3, 'HandleVisibility','off')
plot(trange, [0, 0], ':k', 'HandleVisibility','off')
legend('location', 'southeast')

xlabel('t')
ylabel('y')


%% Nonnegative values
figure(2)
clf
%get data along impulse response trajectory 
xcell = num2cell(xtraj.x, 1);
xtraj.v = pval(xcell{:});
xtraj.Lv = Lpval(xcell{:});
xtraj.cost = (xtraj.x * sys.C').^2;
    

   subplot(3,1,1)
    
    hold on        
    plot(xtraj.t, xtraj.v + obj, 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
    
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')

    hold off
    title('Safety Function along Trajectories')
    xlabel('time')
    ylabel('v(x) - \gamma')
    legend('location', 'northwest')
    
    subplot(3,1,2)
    hold on
    plot(xtraj.t, xtraj.Lv, 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')

    hold off
    title('Change in Safety Function along Trajectories')
    xlabel('time')
    ylabel('L_f v(x)')
    legend('location', 'northwest')
    
    subplot(3,1,3)    
    hold on
    
    p_comp = -xtraj.cost - xtraj.v;
    plot(xtraj.t, p_comp, 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
    
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')
    hold off
    
    
    title('Comparision with Objective Function')
    xlabel('time')
    ylabel('-cost(x) - v(x,t)')
    legend('location', 'northwest')
 

end