function [] = plot_arc(A, B, C, R, anticlockwise)
%PLOT_ARC

   % anticlockwise = 1; clockwise = -1
   % A : Point; B: Candidate, C: Center
   
xCenter = C(1);
yCenter = C(2); 
radius = R;

anticlockwise = anticlockwise*(-1); % Reverse the direction

Ox = [1 0]
a = angle_between(Ox, [B(1) - C(1), B(2) - C(2)])
b = angle_between(Ox, [A(1) - C(1), A(2) - C(2)])

u = [B(1) - C(1), B(2) - C(2)];
v = [A(1) - C(1), A(2) - C(2)];

angle = atan2(v(2), v(1)) - atan2(u(2), u(1));

ang_u_v = (atan2(v(2), v(1)) - atan2(u(2), u(1)))*360/(2*pi)
ang_Ox_u = (atan2(u(2), u(1)) - atan2(Ox(2), Ox(1)))*360/(2*pi)
ang_Ox_v = (atan2(v(2), v(1)) - atan2(Ox(2), Ox(1)))*360/(2*pi)
ang_Ox_u - ang_u_v

count = 1;
theta = [];
if(anticlockwise == -1)
    if(ang_Ox_u > 0)
        if(ang_Ox_u - ang_Ox_v > 0)
            for i=ang_Ox_u:-0.01:(ang_Ox_v )
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
            for i=ang_Ox_u:-0.01:(ang_Ox_v )
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
            for i=ang_Ox_u:0.01:(ang_Ox_v )
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
            for i=ang_Ox_u:0.01:(ang_Ox_v)
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
% if(anticlockwise == -1)
%     a = a + 360
%     b = b + 360
% else
%     a = -a + 360
%     b = -b + 360
% end

% theta = linspace(a, b, 100);

% Define x and y
x = radius * cosd(theta) + xCenter; 
y = radius * sind(theta) + yCenter; 

% [m n] = size(x)
% for i=1:n
%     x(1,i)
%     y(1,i)
%     plot(x(1,i), y(1,i), '.g'); 
%     pause(0.0000000000000000000000000001)
% end
plot(x, y, 'g-','LineWidth', 2); 
end

