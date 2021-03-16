%Code authored by:
% Didier Henrion
% Martine Ganet-Schoeller
% Samir Bennani
%
%
%Run the problem in https://arxiv.org/pdf/1205.2168.pdf
%to see if I get the same solution

I = 27500; % inertia
kp = 2475; kd = 19800; % controller gains
L = 380; % input saturation level
dz1 = 0.2*pi/180; dz2 = 0.05*pi/180; % deadzone levels
thetamax = 50; omegamax = 5; % bounds on initial conditions
epsilon = sqrt(1e-5); % bound on norm of terminal condition
T = 50; % final time
d = input('order of relaxation ='); d = 2*d;
% measures
mpol('x1',2); m1 = meas(x1); % linear regime
mpol('x2',2); m2 = meas(x2); % upper sat
mpol('x3',2); m3 = meas(x3); % lower sat
mpol('x0',2); m0 = meas(x0); % initial
mpol('xT',2); mT = meas(xT); % terminal
% dynamics on normalized time range [0,1]
% saturation input y normalized in [-1,1]
K = -[kp kd]/L;
y1 = K*x1; f1 = T*[x1(2); L*y1/I]; % linear regime
y2 = K*x2; f2 = T*[x2(2); L/I]; % upper sat
y3 = K*x3; f3 = T*[x3(2); -L/I]; % lower set
% test functions for each measure = monomials
g1 = mmon(x1,d); g2 = mmon(x2,d); g3 = mmon(x3,d);
g0 = mmon(x0,d); gT = mmon(xT,d);

% unknown moments of initial measure
y0 = mom(g0);
% unknown moments of terminal measure
yT = mom(gT);
% input LMI moment problem
cost = mom(xT'*xT);
Ay = mom(diff(g1,x1)*f1)+...
mom(diff(g2,x2)*f2)+...
mom(diff(g3,x3)*f3); % dynamics
% trajectory constraints
X = [y1^2<=1; y2>=1; y3<=-1];
% initial constraints
X0 = [x0(1)^2<=thetamax^2, x0(2)^2<=omegamax^2];
% terminal constraints
XT = [xT'*xT<=epsilon^2];
% bounds on trajectory
B = [x1'*x1<=4; x2'*x2<=4; x3'*x3<=4];
% input LMI moment problem
P = msdp(max(cost), ...
mass(m1)+mass(m2)+mass(m3)==1, ...
mass(m0)==1, ...
Ay==yT-y0, ...
X, X0, XT, B);
% solve LMI moment problem
[status,obj] = msol(P);


scale(x1, 4);
scale(x2, 4);
scale(x3, 4);
[status_scale,obj_scale] = msol(P);