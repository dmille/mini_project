clear, close all;

[p,t] = loadmesh('sfterrain.off');

p = p';

n_points = 10;
t_sampling = 1000;
tension = 0;
speed = 500;

waypoints = rand(n_points,2);

waypoints = [waypoints ones(n_points,1) * 0.5];
ordered_pts = shortestpath_through_all_pts(waypoints);

curve = get_curve_data_catmullrom(ordered_pts, t_sampling, tension)';

figure;
pv = subplot(1,2,2);
plane_view(p,t');
hold on;

mmv = subplot(1,2,1);
minimap(p,t');
hold on;
plot3(curve(:,1), curve(:,2), curve(:,3), 'r-');

fly_along_curve(curve, pv, mmv, speed);
%rzview('on')
