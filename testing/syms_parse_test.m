
%symbolic toolbox
syms x y
p = 2*x*y^3 + 3*x^2 -4*x + 2;
p2 = 4*x*(x+3*y + y*(2*x + 2) + 4);

q = 3*x + 4*x^5 - 43;
c = children(p);
c2 = children(p2);
c2s = children(expand(p2));

%conditional expressions
cond = (x <= y + 1);
ccond = children(cond);



%substitute in gloptipoly
mpol('x2')
mpol('y2')
pf = matlabFunction(p);
%ppol = subs(p, [x; y], [x2; y2])
ppol = pf(x2, y2);

condf = matlabFunction(cond);
condpol = condf(x2, y2);