function [curve] = get_curve_data_catmullrom(point, n, tension)
    % duplicate 2 end points
    p_tmp = [point(1,1) point(1,2) point(1,3)];
    point = [p_tmp; point];
    p_tmp = [point(length(point(:,1)),1) point(length(point(:,1)),2) point(length(point(:,1)),3)];
    point = [point; p_tmp];
    
    curve = [];
    for i=1:length(point)-3
        P0 = [point(i,1) point(i,2) point(i,3)];
        P1 = [point(i+1,1) point(i+1,2) point(i+1,3)];
        P2 = [point(i+2,1) point(i+2,2) point(i+2,3)];
        P3 = [point(i+3,1) point(i+3,2) point(i+3,3)];
        
        curve = [curve catmull_rom(P0, P1, P2, P3, n, tension)];
    end
end

