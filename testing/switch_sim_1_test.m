x0 = [1.5; 0];
Tmax = 7;
mu = 1;

rng(50, 'twister')

%[out_sim] = switch_sim(out.dynamics, x0, Tmax, mu, @ode15s);
[out_sim] = switch_sim(out.dynamics, x0, Tmax, mu, @ode45);

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