x0 = [1.5; 0];
Tmax = 7;
mu = 1;

[out_sim] = switch_sim(out.dynamics, x0, Tmax, mu, @ode15s);

plot(out_sim.x(:, 1), out_sim.x(:, 2))