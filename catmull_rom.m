function [curve_data] = catmull_rom(p0, p1, p2, p3, n, tension)
    
%equation P(t) = t_vec*M*point_vec;
    %curve data
    curve_data = [];
    %n: number of interval in 1 spline
    % Parametrization: unit length

    %   dist = norm(p3 - p0);
    %    dt = 1/n*dist;
%   dt = 1/n;
    dt = linspace(0,1,round(n));
    %coeff involving the tension
    k = (1-tension)/2; 
    % matrix after derivation
    M = [0    1    0       0; 
        -k    0    k       0; 
        2*k   k-3  3-2*k  -k;
        -k    2-k  k-2     k];
    
    for i=1:n
        t = dt(i);
        %cubic parameter
        t_vec = [1 t t^2 t^3];
        
        for j=1:3
            point_vec(1,1) = p0(1);
            point_vec(2,1) = p1(1);
            point_vec(3,1) = p2(1);
            point_vec(4,1) = p3(1);
            point_vec(1,2) = p0(2);
            point_vec(2,2) = p1(2);
            point_vec(3,2) = p2(2);
            point_vec(4,2) = p3(2);
            point_vec(1,3) = p0(3);
            point_vec(2,3) = p1(3);
            point_vec(3,3) = p2(3);
            point_vec(4,3) = p3(3);

            temp(j) = t_vec*M*point_vec(:,j);
        end

        for j=1:3
            curve_data(j,i+1) = temp(j);
        end
    end
end

