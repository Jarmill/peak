%I think the reason my scaling code is broken is because I am not properly
%accounting for the coordinate change.

%Time to test it out

mpol('xp', 2, 1);

%support
Xsupp = [];

%dynamics
f1 = [xp(2); -xp(1) + (1/3).* xp(1).^3 - xp(2)];
X1 = xp'*xp <= 4;

%box = [-4, 4; -4, 4];
box = 8;
[~, box_center, box_half] = box_process(2, box);

xp_scale = box_half.*xp + box_center;
xp_inv_scale = (xp - box_center) .* (1./box_half);

f2 = subs((1./box_half) .* f1, xp, xp_scale);

f3 = subs(box_half.*f2, xp, xp_inv_scale);


%test dynamics
C0 = [1.5;0];
C0_scale = eval(xp_inv_scale, xp, C0);
f10 = eval(f1, xp, C0);
f20 = eval(f2, xp, C0_scale);
f2p = f20 .* box_half;

%test support
X2 = subs(X1, xp, xp_scale);
X10 = eval(X1, xp, C0);
X20 = eval(X2, xp, C0);
X2p = eval(X2, xp, C0_scale);