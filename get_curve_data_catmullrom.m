function [curve] = get_curve_data_catmullrom(points, n, tension)
    % duplicate 2 end pointss
    p_tmp = [points(1,1) points(1,2) points(1,3)];
    points = [p_tmp; points];
    p_tmp = [points(length(points(:,1)),1) points(length(points(:,1)),2) points(length(points(:,1)),3)];
    points = [points; p_tmp];
    
    curve = [];
    for i=1:length(points)-3
        P0 = [points(i,1) points(i,2) points(i,3)];
        P1 = [points(i+1,1) points(i+1,2) points(i+1,3)];
        P2 = [points(i+2,1) points(i+2,2) points(i+2,3)];
        P3 = [points(i+3,1) points(i+3,2) points(i+3,3)];
        
        curve = [curve catmull_rom(P0, P1, P2, P3, n, tension)];
    end
end

