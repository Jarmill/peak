%requires running peak_dual_test_alt.m first


%right now restricted to order 3

%aim: recover the function v(x) from the moment program

%vv = monolist([Tc;Xc;Yc],d);
syms tc xc yc;
%vv = monolist([tc; xc; yc], d);
%vv2 = monolist([tc; xc; yc], d+1);

vv = monolist([tc/T; xc; yc], d);
p = dual_rec'*vv + obj;
%p = dual_rec'*vv;

%SOS from yalmip
%p = value(cv)'*vv - value(gamma);

%Moment dual from gloptipoly
%p = dual'*vv + obj;

%export SDP from gloptipoly

%but what is v(x)?

%[xs, ys, infos] = sedumi(model.A, model.b, model.c, model.K);
%zs = model.c - model.A'*ys;
%dual_var = ys(1:length(v));
%dual_var2 = -xs(1:length(vv));
%dual_var = dual_rec; 
%p = dual_var'*vv + obj;
%p = dual_var'*vv;
%p = dual_var'*vv;



%recover the moment matries from dual variables of sedumi export
% s_offset = model.K.f +cumsum( [0, model.K.s.^2]);
% 
% Ms = reshape(zs(s_offset(2) + (1:model.K.s(2)^2)),model.K.s(2),model.K.s(2));
% Mps = reshape(zs(s_offset(3) + (1:model.K.s(3)^2)),model.K.s(3),model.K.s(3));
% 
% %M0 has a redundancy since t==0. therefore the support is reduced.
% 
% %for order 3:
% screen0= [1, 3, 4, 8, 9, 10, 17, 18, 19, 20];
% M0_screen = M0(screen0, screen0);
% ic = 1;
% M0s = reshape(zs(s_offset(ic) + (1:model.K.s(ic)^2)),model.K.s(ic),model.K.s(ic));
% %norm(M0s - M0_screen, 'fro')
% 

% figure(2)
% subplot(1,2,1)
% imagesc(M)
% subplot(1,2,2)
% imagesc(Mzs)


%[pc, pt] = coeffs(p, [tc;xc;yc]);
%pc = double(pc);

%Nt = 200;%
%N = 80;

%[Tc, Xc, Yc] = meshgrid(linspace(0, T, Nt), linspace(-m, m, N), linspace(-m, m, N));

%Vc = subs(p, {tc, xc, yc}, {Tc, Xc, Yc}); 
% 

figure(4)
clf
hold on
plot3(zeros(Ntheta, 1), X0_pts(1, :), X0_pts(2, :), 'k', 'Linewidth', 3)
    
for i = 1:Nsample        
    if i == 1
        plot3(xtraj{i}.t, xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c')
    else
        plot3(xtraj{i}.t, xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'c', 'HandleVisibility','off')
    end        
end  

MD = 60;
fi = fimplicit3(p, [0, T, -m, m, -m, m], 'MeshDensity',MD, ...
    'EdgeColor', 'None','FaceColor', 'k', 'FaceAlpha', 0.3, ...
    'DisplayName', 'Escape Contour');


patch('XData',[0, 0, 1, 1], 'YData',[-m, m, m, -m], 'ZData', gamma_val*[1,1,1,1], 'FaceColor', 'r', 'FaceAlpha', 0.3)

xlim([0, 1])
ylim([-m, m])
zlim([-m, m])
xlabel('scaled time')
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
