function [tout,tk] = chord_length_parametrization(points)
%chord_length_parametrization Summary of this function goes here
%   Detailed explanation goes here

    n = length(points);
    t = zeros(1,n);
    s = 0;

    diff = points - circshift(points,-1); % calculate the differences between x(i) and x(i+1)
    v = vecnorm(diff(1:end-1, :)'); % calculate the euclidean norm of each difference
    tout = [0 cumsum(v)/sum(v)]; % produces normalized cumulative sum array and sets the first element to 0 (special case)
    
    tk = [0 0 0 tout tout(end) tout(end) tout(end)];
end

