clear, close all;

[p,t] = loadmesh('sfterrain.off');

p = p';

n_points = 500;
s = @(x) 0.5 * sin(x * 2 * pi) + 0.5;
l = linspace(0,1, n_points);
f = s(l);
test_curve = [l' f' (ones(n_points,1) * 0.3)]

figure;
pv = subplot(1,2,2);
plane_view(p,t');
hold on;
plot3(test_curve(:,1), test_curve(:,2), test_curve(:,3), '-');

mmv = subplot(1,2,1);
minimap(p,t');
hold on;
plot3(test_curve(:,1), test_curve(:,2), test_curve(:,3), '-');

speed = 100;
fly_along_curve(test_curve, pv, mmv, speed);
%rzview('on')
