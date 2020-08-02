function [x] = trig_sampler(angle_bounds, rect_bounds)
%TRIG_SAMPLER Sample from uniform distribution some angles and rectangular
%             coordinates. angles in angle_bounds, rects in rect_bounds
%output: x = [x_angle; x_rect];

[n_angle, bnd_angle] = size(angle_bounds, 1);
[n_rect, bnd_rect]   = size(rect_bounds,  1);

if bnd_angle == 1
    angle_bounds = [-angle_bounds angle_bounds];
end

if bnd_rect == 1
    rect_bounds = [-rect_bounds rect_bounds];
end



end

