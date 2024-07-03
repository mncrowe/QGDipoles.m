function [F_min, F_max] = F_minmax(F)
% outputs the minimum and maximum values of a field, F

F_min = min(F,[],"all");
F_max = max(F,[],"all");

end

