clear, close all;


%% Begin Global Vars

% Controls the scaling of the terrain mesh in the z axis
zscaling=0.2;

% number of sample points along the curve
t_sampling = 1000;

% Controls something in catmull rom
tension = -0.5;

% Speed is defined as the delay between successive frames (ie
% tsampling points) delay=1/speed
speed = 100;

% Vars to control ribbon width/normal offset
ribbon_width = 0.005;
ribbon_offset = -0.005;

% All waypoints Customize it later
waypoints = rand(10,3);

%% End Global Variables

[p,t] = loadmesh('sfterrain.off', zscaling);
p = p';

n_waypoints = size(waypoints, 1);
%waypoints = [waypoints ones(n_waypoints,1) * 0.5];
ordered_pts = shortestpath_through_all_pts(waypoints);

curve = get_curve_data_catmullrom(ordered_pts, t_sampling, tension)';

figure;
%pv = subplot(1,2,2);
pv = subplot(1,1,1);
plane_view(p,t');
hold on;

plot_path_ribbon(pv, curve, ribbon_width, ribbon_offset);

figure;
%mmv = subplot(1,2,1);
mmv = subplot(1,1,1);
minimap(p,t');
hold on;
plot3(curve(:,1), curve(:,2), curve(:,3), 'r.');

for i = [1:n_waypoints]
    
    plot3(mmv,waypoints(i,1),waypoints(i,2),waypoints(i,3),'g.', ...
                         'MarkerSize', 12);
    [x,y,z] = sphere();
    r = 0.004;
    surf(pv,x*r+waypoints(i,1), y*r+waypoints(i,2), z*r+waypoints(i,3));
    %    plot3(pv,waypoints(i,1),waypoints(i,2),waypoints(i,3),'g.', ...
    %                    'MarkerSize', 40);

end

minimap_cam_height = 1;
fly_along_curve(curve, pv, mmv, speed, minimap_cam_height);
%rzview('on')
