function [path] = dubin_path(point_set, maximum_steering_angle)
point = point_set;
n = length(point(1,:));

figure(1)
clf
plot(point(1,:), point(2,:), '.k');
grid on
axis equal
hold on
plot(point(1,:), point(2,:),'--b');
plot(point(1,:), point(2,:), 'ob', 'MarkerSize', 25);
plot(point(1,:), point(2,:), 'ob', 'MarkerSize', 10);
plot(point(1,:), point(2,:), '.b', 'MarkerSize', 20);
R = maximum_steering_angle;


center = zeros(2,n);
in = zeros(2,n-2);
out = zeros(2,n-2);
out(1,n-2) = point(1,n-1);
out(2,n-2) = point(2,n-1);
direction = zeros(1,n-2);
%% Create closest cirlces
for i=2:n-1
    line = ortho_line_param([point(1,i), point(2,i)], [point(1,i+1), point(2,i+1)]);
    
    if(line(4) == 1)        
        x_temp = point(1,i) + R;
        y_temp = -line(3);
        center1 = [x_temp, y_temp];

        x_temp = point(1,i) - R;
        y_temp = -line(3);
        center2 = [x_temp, y_temp];
        
    elseif(line(4) == 2)        
        x_temp = -line(3);
        y_temp = point(2,i) + R;
        center1 = [x_temp, y_temp];

        x_temp = -line(3);
        y_temp = point(2,i) - R;
        center2 = [x_temp, y_temp];
    else        
        x_temp = sqrt(R*R/((-(line(1))*(-line(1))) + 1)) + point(1,i);
        y_temp = (-line(1))*x_temp + (line(1))*point(1,i) + point(2,i);
        center1 = [x_temp, y_temp];

        x_temp = -sqrt(R*R/((-(line(1))*(-line(1))) + 1)) + point(1,i);
        y_temp = (-line(1))*x_temp + (line(1))*point(1,i) + point(2,i);
        center2 = [x_temp, y_temp];
    end
    
    point_temp = [point(1,i-1), point(2,i-1)];
    if(dist(point_temp, center1) >= dist(point_temp, center2))
        center_x = center2(1);
        center_y = center2(2);
    else
        center_x = center1(1);
        center_y = center1(2);
    end
        
    center(1,i) = center_x;
    center(2,i) = center_y;
end
 
%% Find tangent intersections with cirlces
intersect_store = zeros(2,n);
for i=1:1
    %Tangent intersection
    nominator = R*R*(point(1,i) - center(1,i+1))...
        + R*(point(2,i) - center(2,i+1))*...
        sqrt((point(1,i) - center(1,i+1))*(point(1,i) - center(1,i+1))...
              + (point(2,i) - center(2,i+1))*(point(2,i) - center(2,i+1)) - R*R);
    denominator = (point(1,i) - center(1,i+1))*(point(1,i) - center(1,i+1)) + (point(2,i) - center(2,i+1))*(point(2,i) - center(2,i+1));
    intersect1_x = (nominator)/(denominator) + center(1,i+1);
    
    nominator = R*R*(point(2,i) - center(2,i+1))...
        - R*(point(1,i) - center(1,i+1))*...
        sqrt((point(1,i) - center(1,i+1))*(point(1,i) - center(1,i+1))...
              + (point(2,i) - center(2,i+1))*(point(2,i) - center(2,i+1)) - R*R);
    denominator = (point(1,i) - center(1,i+1))*(point(1,i) - center(1,i+1)) + (point(2,i) - center(2,i+1))*(point(2,i) - center(2,i+1));
    intersect1_y = (nominator)/(denominator) + center(2,i+1);
        
    nominator = R*R*(point(1,i) - center(1,i+1))...
        - R*(point(2,i) - center(2,i+1))*...
        sqrt((point(1,i) - center(1,i+1))*(point(1,i) - center(1,i+1))...
              + (point(2,i) - center(2,i+1))*(point(2,i) - center(2,i+1)) - R*R);
    denominator = (point(1,i) - center(1,i+1))*(point(1,i) - center(1,i+1)) + (point(2,i) - center(2,i+1))*(point(2,i) - center(2,i+1));
    intersect2_x = (nominator)/(denominator) + center(1,i+1);
    
    nominator = R*R*(point(2,i) - center(2,i+1))...
        + R*(point(1,i) - center(1,i+1))*...
        sqrt((point(1,i) - center(1,i+1))*(point(1,i) - center(1,i+1))...
              + (point(2,i) - center(2,i+1))*(point(2,i) - center(2,i+1)) - R*R);
    denominator = (point(1,i) - center(1,i+1))*(point(1,i) - center(1,i+1)) + (point(2,i) - center(2,i+1))*(point(2,i) - center(2,i+1));
    intersect2_y = (nominator)/(denominator) + center(2,i+1);
end

%% Find common tangent of 2 circles
traverse_intersect = zeros(2,n);

for i=2:n-2
    temp1 = (-center(1,i) + center(1,i+1))/2; 
    temp2 = (-center(2,i) + center(2,i+1))/2;
    traverse_intersect(1,i) = center(1,i) + temp1;
    traverse_intersect(2,i) = center(2,i) + temp2;
end

% Establish 4 groups of tangents
direct_g1 = zeros(2,2); direct_g2 = zeros(2,2);
traverse_g1 = zeros(2,2); traverse_g2 = zeros(2,2);

for i=2:n-2      
        %% Find 2 traverse tangents of 2 circles
        a = (center(1,i+1) - center(1,i))^2 - 4*R^2
        b = 2*(center(1,i+1) - center(1,i))*(center(2,i) - center(2,i+1))
        c = (center(2,i) - center(2,i+1))^2 - 4*R^2

        m1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
        m2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
        
        
        bb1 = traverse_intersect(2,i) - m1*traverse_intersect(1,i);
        x = 0:0.01:100;
        y = m1*x + bb1;
%         plot(x,y,'k');
        
        bb2 = traverse_intersect(2,i) - m2*traverse_intersect(1,i);
        x = 0:0.01:100;
        y = m2*x + bb2;
%         plot(x,y,'k');
        
        traverse_line1 = [m1, -1, bb1];
        traverse_line2 = [m2, -1, bb2];
        
        % Find intersections of 2 traverse tangents with 2 circles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a = (traverse_line1(1))^2 + 1
        b = -2*center(1,i) + 2*traverse_line1(1)*(traverse_line1(3) - center(2,i))
        c = (traverse_line1(3) - center(2,i))^2 - R^2 + center(1,i)^2
        
        if (abs((b*b - 4*a*c)) < 0.0000001)
            delta = 0;
        else
            delta = sqrt(b*b - 4*a*c)
        end
        
        x1_traverse = (-b + delta)/(2*a)
        y1_traverse = traverse_line1(1)*x1_traverse + traverse_line1(3)
        
%         plot(x1_traverse, y1_traverse, 'or');

        traverse_g1(1,1) = x1_traverse;
        traverse_g1(2,1) = y1_traverse;
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a = (traverse_line2(1))^2 + 1
        b = -2*center(1,i) + 2*traverse_line2(1)*(traverse_line2(3) - center(2,i))
        c = (traverse_line2(3) - center(2,i))^2 - R^2 + center(1,i)^2
        
        if (abs((b*b - 4*a*c)) < 0.0000001)
            delta = 0;
        else
            delta = sqrt(b*b - 4*a*c)
        end
        
        x2_traverse = (-b + delta)/(2*a)
        y2_traverse = traverse_line2(1)*x2_traverse + traverse_line2(3)
        
        traverse_g2(1,1) = x2_traverse;
        traverse_g2(2,1) = y2_traverse;
%         plot(x2_traverse, y2_traverse, 'or');
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a = (traverse_line1(1))^2 + 1
        b = -2*center(1,i+1) + 2*traverse_line1(1)*(traverse_line1(3) - center(2,i+1))
        c = (traverse_line1(3) - center(2,i+1))^2 - R^2 + center(1,i+1)^2
        
        if (abs((b*b - 4*a*c)) < 0.0000001)
            delta = 0;
        else
            delta = sqrt(b*b - 4*a*c)
        end
        
        x3_traverse = (-b + delta)/(2*a)
        y3_traverse = traverse_line1(1)*x3_traverse + traverse_line1(3)
        
%         plot(x3_traverse, y3_traverse, 'or');
        traverse_g1(1,2) = x3_traverse;
        traverse_g1(2,2) = y3_traverse;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a = (traverse_line2(1))^2 + 1
        b = -2*center(1,i+1) + 2*traverse_line2(1)*(traverse_line2(3) - center(2,i+1))
        c = (traverse_line2(3) - center(2,i+1))^2 - R^2 + center(1,i+1)^2
        
        if (abs((b*b - 4*a*c)) < 0.0000001)
            delta = 0;
        else
            delta = sqrt(b*b - 4*a*c)
        end
        
        x4_traverse = (-b + delta)/(2*a)
        y4_traverse = traverse_line2(1)*x4_traverse + traverse_line2(3)
        
        traverse_g2(1,2) = x4_traverse;
        traverse_g2(2,2) = y4_traverse;
%         plot(x4_traverse, y4_traverse, 'or');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        traverse_point = [x1_traverse x2_traverse x3_traverse x4_traverse; y1_traverse y2_traverse y3_traverse y4_traverse];
        
        %% Find 2 direct tangents of 2 circles
        m = (center(2,i) - center(2,i+1))/(center(1,i) - center(1,i+1));
        bbb1 = -sqrt(m^2+1)*R - m*center(1,i+1) + center(2,i+1);
        bbb2 = sqrt(m^2+1)*R - m*center(1,i+1) + center(2,i+1);
        
        direct_tangent1 = [m, -1, bbb1];
        direct_tangent2 = [m, -1, bbb2];
        
        x = 0:0.01:100;
        y = m*x + bbb1;
%         plot(x,y,'g');
        
        x = 0:0.01:100;
        y = m*x + bbb2;
%         plot(x,y,'g');
        
        % Find interesections of 2 direct tangents with 2 circles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a = (direct_tangent1(1))^2 + 1;
        b = -2*center(1,i) + 2*direct_tangent1(1)*(direct_tangent1(3) - center(2,i));
        c = (direct_tangent1(3) - center(2,i))^2 - R^2 + center(1,i)^2;
        
        if (abs((b*b - 4*a*c)) < 0.0000001)
            delta = 0;
        else
            delta = sqrt(b*b - 4*a*c);
        end
        
        x1_direct = (-b + delta)/(2*a);
        y1_direct = direct_tangent1(1)*x1_direct + direct_tangent1(3);
        direct_g1(1,1) = x1_direct;
        direct_g1(2,1) = y1_direct;
%         plot(x1_direct, y1_direct, 'om');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a = (direct_tangent2(1))^2 + 1
        b = -2*center(1,i) + 2*direct_tangent2(1)*(direct_tangent2(3) - center(2,i))
        c = (direct_tangent2(3) - center(2,i))^2 - R^2 + center(1,i)^2
        
        if (abs((b*b - 4*a*c)) < 0.0000001)
            delta = 0;
        else
            delta = sqrt(b*b - 4*a*c)
        end
        
        x2_direct = (-b + delta)/(2*a)
        y2_direct = direct_tangent2(1)*x2_direct + direct_tangent2(3)
        direct_g2(1,1) = x2_direct;
        direct_g2(2,1) = y2_direct;
%         plot(x2_direct, y2_direct, 'om');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a = (direct_tangent1(1))^2 + 1
        b = -2*center(1,i+1) + 2*direct_tangent1(1)*(direct_tangent1(3) - center(2,i+1))
        c = (direct_tangent1(3) - center(2,i+1))^2 - R^2 + center(1,i+1)^2
        
        if (abs((b*b - 4*a*c)) < 0.0000001)
            delta = 0;
        else
            delta = sqrt(b*b - 4*a*c)
        end
        
        x3_direct = (-b + delta)/(2*a)
        y3_direct = direct_tangent1(1)*x3_direct + direct_tangent1(3)
        direct_g1(1,2) = x3_direct;
        direct_g1(2,2) = y3_direct;
%         plot(x3_direct, y3_direct, 'om');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a = (direct_tangent2(1))^2 + 1
        b = -2*center(1,i+1) + 2*direct_tangent2(1)*(direct_tangent2(3) - center(2,i+1))
        c = (direct_tangent2(3) - center(2,i+1))^2 - R^2 + center(1,i+1)^2
        
        if (abs((b*b - 4*a*c)) < 0.0000001)
            delta = 0;
        else
            delta = sqrt(b*b - 4*a*c)
        end
        
        x4_direct = (-b + delta)/(2*a)
        y4_direct = direct_tangent1(1)*x4_direct + direct_tangent2(3)
        direct_g2(1,2) = x4_direct;
        direct_g2(2,2) = y4_direct;
%         plot(x4_direct, y4_direct, 'om');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        direct_point = [x1_direct x2_direct x3_direct x4_direct; y1_direct y2_direct y3_direct y4_direct];
        
        %% Find optimal tangent path
        
%         plot(traverse_g1(1,:), traverse_g1(2,:),'*r')
%         plot(traverse_g2(1,:), traverse_g2(2,:),'*g')
        
%         plot(direct_g1(1,:), direct_g1(2,:),'*m')
%         plot(direct_g2(1,:), direct_g2(2,:),'*y')
        
      
        sign_circle1 = zeros(1,4);   % Structure = (traverse1 traverse2 direct1 direct2)
        sign_circle2 = zeros(1,4);   % Structure = (traverse1 traverse2 direct1 direct2)
        
        Tvec_traverse_g1 = [traverse_g1(1,2) - traverse_g1(1,1), traverse_g1(2,2) - traverse_g1(2,1)]
        Tvec_direct_g1 = [direct_g1(1,2) - direct_g1(1,1), direct_g1(2,2) - direct_g1(2,1)];
        Tvec_traverse_g2 = [traverse_g2(1,2) - traverse_g2(1,1), traverse_g2(2,2) - traverse_g2(2,1)];
        Tvec_direct_g2 = [direct_g2(1,2) - direct_g2(1,1), direct_g2(2,2) - direct_g2(2,1)];
        
        %%%%%Transform target vector to tangent points at circle (i+1)%%%%%
        %%% traverse tangents
        
        target_vec = [point(1,i+2) - point(1,i+1), point(2,i+2) - point(2,i+1)];
        u = target_vec;
        v = -[center(1,i+1) - point(1,i+1), center(2,i+1) - point(2,i+1)]
        KK = [u(1) u(2); v(1) v(2)]
        sign_target_vec = det(KK)/abs(det(KK));
        
        %rotate to match traslate vector in group 1
        u = Tvec_traverse_g1;
        v = -[center(1,i+1) - traverse_g1(1,2), center(2,i+1) - traverse_g1(2,2)]
        KK = [u(1) u(2); v(1) v(2)];
        sign = det(KK)/abs(det(KK))
        sign = sign*sign_target_vec
        sign_circle2(1) = sign;
        theta = angle_between(target_vec, sign*Tvec_traverse_g1)
        Ro = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]
        target_new21_traverse = Ro*target_vec'
        angle_between(target_new21_traverse', sign*Tvec_traverse_g1)
                
        %rotate to match traslate vector in group 2
        u = Tvec_traverse_g2;
        v = -[center(1,i+1) - traverse_g2(1,2), center(2,i+1) - traverse_g2(2,2)];
        KK = [u(1) u(2); v(1) v(2)];
        sign = det(KK)/abs(det(KK));
        sign = sign*sign_target_vec
        sign_circle2(2) = sign;
        theta = angle_between(target_vec, sign*Tvec_traverse_g2)
        Ro = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]
        target_new22_traverse = Ro*target_vec'
        angle_between(target_new22_traverse', sign*Tvec_traverse_g2)
        
        %%% direct tangent
        %rotate to match traslate vector in group 1
        u = Tvec_direct_g1;
        v = -[center(1,i+1) - direct_g1(1,2), center(2,i+1) - direct_g1(2,2)];
        KK = [u(1) u(2); v(1) v(2)];
        sign = det(KK)/abs(det(KK));
        sign = sign*sign_target_vec
        sign_circle2(3) = sign;
        theta = angle_between(target_vec, sign*Tvec_direct_g1)
        Ro = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]
        target_new21_direct = Ro*target_vec'
        angle_between(target_new21_direct', sign*Tvec_direct_g1)
                
        %rotate to match traslate vector in group 2
        u = Tvec_direct_g2;
        v = -[center(1,i+1) - direct_g2(1,2), center(2,i+1) - direct_g2(2,2)];
        KK = [u(1) u(2); v(1) v(2)];
        sign = det(KK)/abs(det(KK));
        sign = sign*sign_target_vec
        sign_circle2(4) = sign;
        theta = angle_between(target_vec, sign*Tvec_direct_g2)
        Ro = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]
        target_new22_direct = Ro*target_vec'
        angle_between(target_new22_direct', sign*Tvec_direct_g2)

       
        
        %%%%% Translate target vectors to tangent at circle (i) %%%%% 
        %%% traverse tangents
        target_vec = [point(1,i+1) - point(1,i), point(2,i+1) - point(2,i)];
        u = target_vec;
        v = -[center(1,i) - point(1,i), center(2,i) - point(2,i)]
        KK = [u(1) u(2); v(1) v(2)]
        sign_target_vec = det(KK)/abs(det(KK));
        
        %rotate to match traslate vector in group 1
        u = Tvec_traverse_g1;
        v = -[center(1,i) - traverse_g1(1,1), center(2,i) - traverse_g1(2,1)];
        KK = [u(1) u(2); v(1) v(2)];
        sign = det(KK)/abs(det(KK));
        sign = sign*sign_target_vec
        sign_circle1(1) = sign;
        theta = angle_between(target_vec, sign*Tvec_traverse_g1)
        Ro = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]
        target_new11_traverse = Ro*target_vec'
        angle_between(target_new11_traverse', sign*Tvec_traverse_g1)
        
        %rotate to match traslate vector in group 2
        u = Tvec_traverse_g2;
        v = -[center(1,i) - traverse_g2(1,1), center(2,i) - traverse_g2(2,1)];
        KK = [u(1) u(2); v(1) v(2)];
        sign = det(KK)/abs(det(KK));
        sign = sign*sign_target_vec
        sign_circle1(2) = sign;
        theta = angle_between(target_vec, sign*Tvec_traverse_g2)
        Ro = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]
        target_new12_traverse = Ro*target_vec'
        angle_between(target_new12_traverse', sign*Tvec_traverse_g2)
        
        %%% direct tangent
        %rotate to match traslate vector in group 1
        u = Tvec_direct_g1;
        v = -[center(1,i) - direct_g1(1,1), center(2,i) - direct_g1(2,1)];
        KK = [u(1) u(2); v(1) v(2)];
        sign = det(KK)/abs(det(KK));
        sign = sign*sign_target_vec;
        sign_circle1(3) = sign;
        theta = angle_between(target_vec, sign*Tvec_direct_g1);
        Ro = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
        target_new11_direct = Ro*target_vec';
        angle_between(target_new11_direct', sign*Tvec_direct_g1)
                
        %rotate to match traslate vector in group 2
        u = Tvec_direct_g2;
        v = -[center(1,i) - direct_g2(1,1), center(2,i) - direct_g2(2,1)];
        KK = [u(1) u(2); v(1) v(2)];
        sign = det(KK)/abs(det(KK));
        sign = sign*sign_target_vec
        sign_circle1(4) = sign;
        theta = angle_between(target_vec, sign*Tvec_direct_g2)
        Ro = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]
        target_new12_direct = Ro*target_vec'
        angle_between(target_new12_direct', sign*Tvec_direct_g2)
        
        
        %%%%%Find optimal path%%%%%     
        sign_circle1
        sign_circle2
        sign_combine = sign_circle1.*sign_circle2
        
        cell = {traverse_g1; traverse_g2; direct_g1; direct_g2};
        candidate = zeros(2,4);
        ind = 1;
        marker = zeros(1,4);
        for j=1:4
            if(sign_combine(j) == 1)
                candidate(1,ind) = cell{j}(1,1)
                candidate(2,ind) = cell{j}(2,1)
                candidate(1,ind+1) = cell{j}(1,2)
                candidate(2,ind+1) = cell{j}(2,2)
                ind = ind + 2
            end
        end
        
        %%%%%Find best fit vector pair%%%%%%
%         plot(candidate(1,:), candidate(2,:), '.k', 'MarkerSize', 10)
        
        target_vec = [point(1,i+1) - point(1,i), point(2,i+1) - point(2,i)];
        u = target_vec;
        v = -[center(1,i) - point(1,i), center(2,i) - point(2,i)]
        KK = [u(1) u(2); v(1) v(2)];
        sign_target_vec = det(KK)/abs(det(KK))
        
        
        tangent1_vec = [candidate(1,2) - candidate(1,1), candidate(2,2) - candidate(2,1)];
        u = tangent1_vec;
        v = -[center(1,i) - candidate(1,2), center(2,i) - candidate(2,2)];
        KK = [u(1) u(2); v(1) v(2)];
        sign1 = det(KK)/abs(det(KK));
        
        tangent2_vec = [candidate(1,4) - candidate(1,3), candidate(2,4) - candidate(2,3)];
        u = tangent2_vec;
        v = -[center(1,i) - candidate(1,4), center(2,i) - candidate(2,4)];
        KK = [u(1) u(2); v(1) v(2)];
        sign2 = det(KK)/abs(det(KK));
        
        if(sign1*sign_target_vec == 1)
            plot(candidate(1,1:2), candidate(2,1:2), '-k','LineWidth', 2)
            
            in(1, i) = candidate(1,2);
            in(2, i) = candidate(2,2);
            
            out(1, i-1) = candidate(1,1);
            out(2, i-1) = candidate(2,1);
            
            temp_vec1 = [candidate(1,2) - candidate(1,1), candidate(2,2) - candidate(2,1)];
            u = temp_vec1;
            v = -[center(1,i) - candidate(1,2), center(2,i) - candidate(2,2)];
            KK = [u(1) u(2); v(1) v(2)];
            sign1 = det(KK)/abs(det(KK));
            direction(i-1) = sign1;
        else
            plot(candidate(1,3:4), candidate(2,3:4), '-k','LineWidth', 2)
            
            in(1, i) = candidate(1,4);
            in(2, i) = candidate(2,4);
            
            out(1, i-1) = candidate(1,3);
            out(2, i-1) = candidate(2,3);
            
            temp_vec2 = [candidate(1,4) - candidate(1,3), candidate(2,4) - candidate(2,3)];
            u = temp_vec2 ;
            v = -[center(1,i) - candidate(1,4), center(2,i) - candidate(2,4)];
            KK = [u(1) u(2); v(1) v(2)];
            sign2 = det(KK)/abs(det(KK));
            direction(i-1) = sign2;
        end
end 

% % % % % % Process 1st point
    temp_vec = [in(1,2) - out(1,1),in(2,2) - out(2,1)]
    u = temp_vec ;
    v = -[center(1,2) - out(1,1),center(2,2) - out(2,1)]
    KK = [u(1) u(2); v(1) v(2)]
    sign = det(KK)/abs(det(KK))

    d1 = dist_evaluate_at_p1([out(1,1) out(2,1)], [intersect1_x intersect1_y], [center(1,2) center(2,2)], R, sign)
    d2 = dist_evaluate_at_p1([out(1,1) out(2,1)], [intersect2_x intersect2_y], [center(1,2) center(2,2)], R, sign)

    if(d1 < d2)
        in(1,1) = intersect1_x
        in(2,1) = intersect1_y;
    else
        in(1,1) = intersect2_x
        in(2,1) = intersect2_y;
    end
    
plot([point(1,1) in(1,1)], [point(2,1) in(2,1)],'-k', 'LineWidth', 2)

% Create plot/data
temp_vec2 = [in(1,n-2) - out(1,n-3), in(2,n-2) - out(2,n-3)];
u = temp_vec2 ;
v = -[center(1,n-1) - in(1,n-2), center(2,n-1) - in(2,n-2)];
KK = [u(1) u(2); v(1) v(2)];
sign2 = det(KK)/abs(det(KK))
direction(n-2) = sign2;
 
plot(point(1,n-1:n), point(2,n-1:n),'-k', 'LineWidth', 2)

path = zeros(2,(n-2)*2+2);
path(1,1) = point(1,1);
path(2,1) = point(2,1);

path(1,(n-2)*2+2) = point(1,n);
path(2,(n-2)*2+2) = point(2,n);
count = 1;

%Pattern  LINE---CURVE---LINE---CURVE---...---LINE
for i=1:n-2
    count = count + 1;
    path(1,count) = in(1,i);
    path(2,count) = in(2,i);
    count = count + 1;
    path(1,count) = out(1,i);
    path(2,count) = out(2,i);
end

plot(point(1,1), point(2,1), 'om', 'MarkerSize', 35);
plot(point(1,1), point(2,1), 'om', 'MarkerSize', 20);
plot(point(1,1), point(2,1), '.m', 'MarkerSize', 30);

plot(point(1,n), point(2,n), 'ok', 'MarkerSize', 35);
plot(point(1,n), point(2,n), 'ok', 'MarkerSize', 20);
plot(point(1,n), point(2,n), '.k', 'MarkerSize', 30);


% GETTING CURVE DATA
curve_data = {};
for i=1:n-2
%      plot_arc([out(1,i) out(2,i)], [in(1,i) in(2,i)],[center(1,i+1) center(2,i+1)], R, direction(i))
    flag = 0;
    if(i < n-2)
        d1 = dist([point(1,i+1) point(2,i+1)], [in(1,i+1) in(2,i+1)])
        d2 = dist([out(1,i) out(2,i)], [in(1,i+1) in(2,i+1)])
        if(d1 < d2)
            flag = 1
        end
    end
    curve_data = [curve_data get_curve_data([out(1,i) out(2,i)], [in(1,i) in(2,i)],[center(1,i+1) center(2,i+1)], R, direction(i), flag)];
end

% plot(curve_data{:}(1,:), curve_data{:}(2,:), '*r')

%GETTING LINE DATA
line_data = {};

x0 = point(1,1);
y0 = point(2,1);
m = -[point(1,1) - in(1,1), point(2,1) - in(2,1)];
line_dat_tmp = [];
k = 0;
for i=0:1/36000:1
    x = x0 + i*m(1);
    y = y0 + i*m(2);
    k = k+1;
    line_dat_tmp(1,k) = x;
    line_dat_tmp(2,k) = y;
end
line_data = [line_data line_dat_tmp];

for j=1:length(in(1,:))-1
    x0 = out(1,j);
    y0 = out(2,j);
    m = [in(1,j+1) - x0, in(2,j+1) - y0];
    line_dat_tmp = [];
    k = 0;
    for i=0:1/36000:1
        x = x0 + i*m(1);
        y = y0 + i*m(2);
        k = k+1;
        line_dat_tmp(1,k) = x;
        line_dat_tmp(2,k) = y;
    end
    line_data = [line_data line_dat_tmp];
end

x0 = out(1,length(out(1,:)));
y0 = out(2,length(out(1,:)));
m = [point(1,n) - x0, point(2,n) - y0];
line_dat_tmp = [];
k = 0;
for i=0:1/36000:1
    x = x0 + i*m(1);
    y = y0 + i*m(2);
    k = k+1;
    line_dat_tmp(1,k) = x;
    line_dat_tmp(2,k) = y;
end

line_data = [line_data line_dat_tmp];

%FINAL PATH DATA
final_data = [];
for i=1:length(curve_data)
    final_data = [final_data line_data{i}];
    final_data = [final_data curve_data{i}];
end

final_data = [final_data line_data{length(line_data)}];

plot(final_data(1,:), final_data(2,:), '.g');

path  = final_data;


test_bspline_fitting;

end