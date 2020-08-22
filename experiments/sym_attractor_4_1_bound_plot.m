peak_ref_fin = [2.062962, 1.918262, 1.901441, 1.901409, 1.901409, 1.901409, 1.901409, 1.901409]';

peak_ref_inf = [2.194343, 1.942396, 1.931330, 1.916228, 1.903525, 1.903448, 1.903185, 1.903181]';

peak_struct = load('experiments/sym_attractor_experiment.mat', 'degree_bound', 'peak_fin_horizon', 'peak_inf_horizon');

figure(1)
clf
hold on
plot(2*degree_bound, peak_ref_fin - peak_struct.peak_fin_horizon)
plot(2*degree_bound, peak_ref_inf - peak_struct.peak_inf_horizon)
set(gca, 'Yscale', 'log')
xlabel('degree of v')
ylabel('Difference in Upper Bounds')