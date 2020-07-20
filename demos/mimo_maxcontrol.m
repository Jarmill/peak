rng(10)

nx = 4;
ny = 2;
nu = 3;

sys = rss(nx, ny, nu);
%X0 = [eye(nx) -eye(nx)];
X0 = [eye(nx) ones(nx, 1)*0.6];
%X0 = [1 0 2/3; 0 1 2/3];

Q = norm(sys.C'*sys.C, 'fro') * eye(nx);
R = 1;

[K, S, e] = lqr(sys, Q, R);

order = 2;
[peak_val, out] = peak_lin_control(sys.A, sys.B, K, X0, order);

