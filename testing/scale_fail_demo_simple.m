% mpol x 1
% mpol y 1
% 
% g0 = -2*x + y;
% 
% K = [x^4 + y^4  - 3*x - 4* y <= 10, x + y >= 1];
% 
% P = msdp(min(g0), K);
% 
% [status, obj] = msol(P);
% 
% scale(x, 5), scale(y, 5)
% K_scale = [x^4 + y^4  - 3*x - 4* y <= 10, x + y >= 1];
% P_scale= msdp(min(g0), K_scale);
% 
% [status_scale, obj_scale] = msol(P_scale);

mpol x
g0 = x^2 - 2*x;
h0 = x^2 + 2*x + 1;
P = msdp(min(g0), mom(h0) == 1);
%scale(x, 100)
[status, obj] = msol(P);