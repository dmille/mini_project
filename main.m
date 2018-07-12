clear, close all;

[p,t] = loadmesh('sfterrain.off');

p = p';

n_points = 10;
t_sampling = 1000;
tension = -0.5;
speed = 100;

waypoints = rand(n_points,3);

%waypoints = [waypoints ones(n_points,1) * 0.5];
ordered_pts = shortestpath_through_all_pts(waypoints);

curve = get_curve_data_catmullrom(ordered_pts, t_sampling, tension)';

figure;
pv = subplot(1,2,2);
plane_view(p,t');
hold on;

%plot_path_ribbon(pv, curve);

plot3(curve(:,1), curve(:,2), curve(:,3)-0.005, ':r');

mmv = subplot(1,2,1);
minimap(p,t');
hold on;
plot3(curve(:,1), curve(:,2), curve(:,3), 'r.');
for i = [1:n_points]
    
    plot3(mmv,waypoints(i,1),waypoints(i,2),waypoints(i,3),'g.', ...
                         'MarkerSize', 12);
    [x,y,z] = sphere();
    r = 0.004;
    surf(pv,x*r+waypoints(i,1), y*r+waypoints(i,2), z*r+waypoints(i,3));
    %    plot3(pv,waypoints(i,1),waypoints(i,2),waypoints(i,3),'g.', ...
    %                    'MarkerSize', 40);

end

fly_along_curve(curve, pv, mmv, speed);
%rzview('on')
