% Estimate moments of a trajectory solving a delay differential equation
% x'(t) = -K x(t - tau) where the state history is constant at xh0
%
% Right now there is a disagreement in the moments, refer to the variable 
% 'm_comp'. Empirical moments from a trajectory are on the left, and
% moments from the LMI (t,x0 marginal of joint occupation measure) are on 
% the right. The (t) marginal is the same, but moments with x0 are not
% aligned. Disagreement may be due to issues in theory or implementation.
% 
% Author: Jared Miller
%         Feb 2, 2021.

%% parameters
PLOT = 1;
T = 1;      %time horizon
tau = 0.1;  %delay x(t - tau)
xh0 = -1;   %constant history x(t) = xh0 for times [-tau, 0]
K = 10;     %gain in dynamics 
order = 4;      %relaxation order

%% plot the trajectory
options = ddeset('AbsTol', 1e-11, 'RelTol', 1e-8);
%ddesd(dynamics, delay, history, time range, options)
sol = ddesd(@(t,y,z) -K*z, @(t,y) t-tau,@(t) xh0,[0,T], options);

if PLOT
    figure(1)
    clf
    hold on
    plot([-tau sol.x], [-1 sol.y])
    plot([-tau, T], [0, 0], ':k')
    xlim([-tau, T])
    hold off
    xlabel('t')
    ylabel('x(t)')
    title(['$\dot{x}(t) = -', num2str(K), 'x(t-', num2str(tau),')$'], 'interpreter', 'latex', 'fontsize', 16)
end

%% Set up variables and measures
d = 2*order;

mset clear
mpol('t','x0','x1');mu  = meas(t, x0, x1);  %joint occupation measure
mpol('tT', 'xT');   muT = meas(tT, xT);     %final measure
mpol('tnn', 'xnn'); nun = meas(tnn, xnn);   %component -1
mpol('tnz', 'xnz'); nuz = meas(tnz, xnz);   %component  0
mpol('tnp', 'xnp'); nup = meas(tnp, xnp);   %component  1

%support constraints
supp_con = [tT == T; 
            tnz*((T - tau) - tnz) >= 0; 
            (tnp - (T - tau)) * (T - tnp) >= 0;
            t*(T - t) >= 0;            
            x0^2 <= xh0^2;
            x1^2 <= xh0^2;
            xT^2 <= xh0^2;
            xnz^2 <= xh0^2;
            xnp^2 <= xh0^2];
%                     xnn == xh0; % redundant
%             (tnn + tau) * tnn >= 0; %redundant
        
%% monomials
yT = mmon([tT, xT], d);
        
%% reference variables and measures
dv = genPowGlopti(2,d);
y0 = prod([0, xh0].^dv, 2); %initial measure
yh = history_mom(tau, xh0, dv);


%% Affine relation
%Liouville Equation
v0  = mmon([t; x0], d);
f = -5*x1; %dynamics
Ay = mom(diff(v0, t) + diff(v0, x0)*f); 
Liou = Ay + y0 - mom(yT);

%(t, x0) marginal
%test functions on components (t)^alpha x^beta
phi0 = mom(v0) - mom(mmon([tnz; xnz], d)) - mom(mmon([tnp; xnp], d));

%(t, x1) marginal (with shifting)
%test functions on components (t+tau)^alpha x^beta
% joint occupation measure
v1  = mmon([t; x1], d);

%[0, T-tau] component
mnz = mmon([tnz; xnz], d);
mnz_shift = subs(mnz, tnz, tnz + tau);  

%[-tau, 0] component]
mnn = mmon([tnn; xnn], d);
mnn_shift = subs(mnn, tnn, tnn + tau);

phi1 = mom(v1) - mom(mnz_shift) - mom(mnn_shift);



%% moment constraints
mom_con = [Liou == 0;        %Liouville
           phi0 == 0;        %x(t)
           phi1 == 0;        %x(t - tau)
           mom(mnn) == yh];  %history given

%% Solve problem
% objective = max(mass(mnn)); 
objective = max(mass(nuz)); 

P = msdp(objective, ...
    mom_con, supp_con);

% solve LMIP moment problem
[status, obj] = msol(P);


%% Analyze moments

%estimation from trajectory
m_traj = monom_int(sol.x, sol.y, dv);
%(t, x0) marginal of joint occupation measure mu
m_mom  = double(mom(v0));
m_comp = [m_traj m_mom];
norm_diff = norm(m_traj-m_mom);

disp(['Norm Gap between empirical and LMI moments = ', num2str(norm_diff)])

function m = history_mom(tau, x0, dv)
    %moments of (t: lebesgue of [-tau, 0]) times (x: delta at x0)
    alpha = dv(:, 1);
    beta  = dv(:, 2);
    
    tmom = -(-tau).^(alpha+1)./(alpha+1);
    xmom = x0.^beta;
    
    m = xmom .* tmom;
end

function em = monom_int(t, x, dv)
    %numerically integrate t^alpha x(t)^beta in span [0, T]
    %em: empirical moments
    %there is probably a more efficient way to implement this
    n_monom = size(dv, 1);
    em = zeros(n_monom, 1);
    for i = 1:n_monom        
        v_curr = (t.^dv(i, 1)) .* (x.^dv(i, 2));
        em(i) = trapz(t, v_curr);
    end
end