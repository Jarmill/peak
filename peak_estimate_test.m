%test out arguments of peak_estimate routine
%clear
%mset clear


%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 2, 1);

%support
Xsupp = [];

%dynamics
f1 = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
f2 = [-x(1); -x(2)];
f3 = [-0.1*x(1)-x(2); x(1) - 0.1*x(2)];
X1 = [];
X2 = [];
X3 = [];
%f = {f1, f2, f3};
%X = {X1, X2, X3};

f = {f1, f3};
X = {X1, X3};
% 
% f = f1;
% X = X1;


%dynamics = struct('f', f, 'X', X);
%dynamics = struct('f', {f1, f2}, 'X', {X, X});


%initial set
C0 = [1.5; 0];
R0 = 0.4;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%objective to maximize
objective = -x(2);
%objective = -x(2) - x(1);
%
p_opt = peak_options;
p_opt.var.x = x;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

%p_opt.Tmax = 10;
p_opt.state_init = X0;
p_opt.state_supp = Xsupp;
p_opt.R = 5;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 6;
out = peak_estimate(p_opt, order);



%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;
Tmax_sim = 15;
mu = 0.5;
[out_sim] = switch_sim(out.dynamics, x0, Tmax_sim, mu, @ode45);

figure(1)
clf
hold on
plot(out_sim.x(:, 1), out_sim.x(:, 2))
scatter(x0(1), x0(2), [], 200, 'ko')

%implicit plot
syms tc xc yc;
if out.dynamics.time_indep
    vv = conj(monolistYalToGlop([xc; yc], 2*order));
end
vsym = out.func.dual_rec'*vv;

fi2 = fimplicit(vsym  + out.peak_val, [-0.5, 1.5, -2, 2], ...
            ':k', 'DisplayName', 'Peak Contour', 'LineWidth', 3, 'MeshDensity', 150);    


figure(2)
clf

nn = out_sim.nonneg;
Nnn = size(out_sim.nonneg, 1);
%tiledlayout(Nnn + 2, 1);
for i = 1:Nnn
    subplot(Nnn, 1, i)
    hold on
    plot(out_sim.t, out_sim.nonneg(i, :));
    plot(out_sim.t([1, end]), [0, 0], ':k')
    hold off
    %nexttile;
end