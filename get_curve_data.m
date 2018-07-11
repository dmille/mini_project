function [curve] = get_curve_data(A, B, C, R, anticlockwise, flag)
%PLOT_ARC

   % anticlockwise = 1; clockwise = -1 (original CCW = -1; CW = 1)
   % A : Point; B: Candidate, C: Center
   
xCenter = C(1);
yCenter = C(2); 
radius = R;

anticlockwise = anticlockwise*(-1); % Reverse the direction

Ox = [1 0];

u = [B(1) - C(1), B(2) - C(2)];
v = [A(1) - C(1), A(2) - C(2)];

% ang_u_v = (atan2(v(2), v(1)) - atan2(u(2), u(1)))*360/(2*pi);
ang_Ox_u = (atan2(u(2), u(1)) - atan2(Ox(2), Ox(1)))*360/(2*pi);
ang_Ox_v = (atan2(v(2), v(1)) - atan2(Ox(2), Ox(1)))*360/(2*pi);

count = 1;
theta = [];
if(anticlockwise == -1)
    if(ang_Ox_u > 0)
        if(ang_Ox_u - ang_Ox_v > 0)
            for i=ang_Ox_u:-0.01:(ang_Ox_v - flag*360)
                theta(count) = i;
                count = count + 1;
            end
        else
            for i=ang_Ox_u:-0.01:-(-ang_Ox_v + 360)
                theta(count) = i;
                count = count + 1;
            end
        end
    else
       if(ang_Ox_u - ang_Ox_v > 0)
            for i=ang_Ox_u:-0.01:(ang_Ox_v -flag*360)
                theta(count) = i;
                count = count + 1;
            end
        else
            for i=ang_Ox_u:-0.01:-(-ang_Ox_v + 360)
                theta(count) = i;
                count = count + 1;
            end
        end
    end
else
    if(ang_Ox_u > 0)
        if(ang_Ox_u - ang_Ox_v < 0)
            for i=ang_Ox_u:0.01:(ang_Ox_v + flag*360)
                theta(count) = i;
                count = count + 1;
            end
        else
            for i=ang_Ox_u:0.01:(ang_Ox_v + 360)
                theta(count) = i;
                count = count + 1;
            end
        end
    else
       if(ang_Ox_u - ang_Ox_v < 0)
            for i=ang_Ox_u:0.01:(ang_Ox_v + flag*360)
                theta(count) = i;
                count = count + 1;
            end
        else
            for i=ang_Ox_u:0.01:(ang_Ox_v + 360)
                theta(count) = i;
                count = count + 1;
            end
        end
    end
end

% Define x and y
x = radius * cosd(theta) + xCenter; 
y = radius * sind(theta) + yCenter; 

[m,n] = size(x);
curve = zeros(2, n);
for i=1:n
    curve(1,i) = x(i);
    curve(2,i) = y(i);
end
% plot(x, y, 'g-','LineWidth', 2); 
end

