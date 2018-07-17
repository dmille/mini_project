function [tout] = reparam_curve(curve,waypoints)
%reparam_curve 
%   after generating first curve, use this to get new tknot values for make
%   better time intervals. 
    n = length(curve);
    t = zeros(1,n);
    s = 0;

    diff = curve - circshift(curve,-1); % calculate the differences between x(i) and x(i+1)
    v = vecnorm(diff(1:end-1, :)'); % calculate the euclidean norm of each difference
    tout = [0 cumsum(v)/sum(v)]; % produces normalized cumulative sum array and sets the first element to 0 (special case)
    
    t_waypoints = find(ismember(curve,waypoints,'rows') == 1);
    diff = t_waypoints - circshift(t_waypoints,-1);
    diff1 = find(diff == -1);
    t_waypoints(diff1) = [];
    
    tout = tout(t_waypoints);
    tout = circshift(tout,-1) - tout; % calculates the difference between x(i) and x(i+1)
    tout = tout(1:end-1); % cuts off end
    
end

