function [new_lim] = stretch(lim,stretch)
%STRETCH Stretch the existing bound in 'lim' by a factor 'stretch'
%[-1, 1] -> [-2, 2] with a stretch of 2
%[0, 2] -> [-1, 3] with a stretch of 2

new_lim = sum(lim)/2 + diff(lim/2)*[-stretch, stretch];
end

