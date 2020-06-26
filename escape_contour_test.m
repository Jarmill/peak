%requires running peak_dual_test_alt.m first


%vv = monolist([Tc;Xc;Yc],d);
syms tc xc yc;
vv = monolist([tc; xc; yc], d);

%SOS from yalmip
%p = value(cv)'*vv - value(gamma);

%Moment dual from gloptipoly
%p = dual'*vv + obj;

%export SDP from gloptipoly
[xs, ys, infos] = sedumi(model.A, model.b, model.c, model.K);
dual_var = xs(1:length(v));
p = dual_var'*vv - obj;
%[pc, pt] = coeffs(p, [tc;xc;yc]);
%pc = double(pc);

%Nt = 200;%
%N = 80;

%[Tc, Xc, Yc] = meshgrid(linspace(0, T, Nt), linspace(-m, m, N), linspace(-m, m, N));

%Vc = subs(p, {tc, xc, yc}, {Tc, Xc, Yc}); 

MD = 50;
fi = fimplicit3(p, [0, T, -m, m, -m, m], 'MeshDensity',MD, ...
    'EdgeColor', 'None','FaceColor', 'k', 'FaceAlpha', 0.3, ...
    'DisplayName', 'Escape Contour');


xlim([0, T])
ylim([-m, m])
zlim([-m, m])
xlabel('time')
ylabel('x_1')
zlabel('x_2')
pbaspect([2.25 1 1])



%xb = R;
%Nt = 400;
%N = 100;
%[Tc, Xc, Yc] = meshgrid(linspace(0, T, Nt), linspace(-xb, xb, N), linspace(-xb, xb, N));
%V = eval(p);
%contour(X,Y,Z, [1 1], '-b', 'linewidth',2); hold on
%mesh(X, Y, Z)

%isosurface(Tc,Xc,Yc,V,0)
