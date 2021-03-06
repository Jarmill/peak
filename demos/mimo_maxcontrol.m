rng(10)

nx = 4;
ny = 2;
nu = 2;

sys = rss(nx, ny, nu);
X0 = [eye(nx) ones(nx, 1)*0.7];
%X0 = [1 0 2/3; 0 1 2/3];

nx0 = size(X0, 2);
n_problem = nx0 * nu;

Q = norm(sys.C'*sys.C, 'fro') * eye(nx);
R = 1;

[K, S, e] = lqr(sys, Q, R);

order = 2;
Tmax = 10;
%maximum control effort
[peak_val, out] = peak_lin_sys(sys.A - sys.B * K, K, X0, order, 0,  1e-3, Tmax);

%maximum state absolute value
%[peak_val_state, out_state] = peak_lin_sys(sys.A - sys.B * K, eye(nx), X0, order);