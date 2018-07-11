function [angle] = angle_between(u,v)
%ANGLE_BETWEEN Summary of this function goes here
%   Detailed explanation goes here
    if (length(u) < 3)
        u = [u 0];
    end
    if (length(v) < 3)
        v = [v 0];
    end
    
   angle = atan2d(v(2), v(1)) - atan2d(u(2), u(1)) ;
end

