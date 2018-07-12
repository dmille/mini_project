p0 = [0 0 0]
p1 = [1 1 1]
p2 = [2 2 1]
p3 = [3 3 0]
tension = 0.3
n = 100
p=[p0; p1; p2; p3]
te = []

te = get_curve_data_catmullrom(p,n, tension);
figure(1)
clf
grid on
axis equal
hold on
plot3(p0(1),p0(2),p0(3),'ob')
plot3(p1(1),p1(2),p1(3),'ob')
plot3(p2(1),p2(2),p2(3),'ob')
plot3(p3(1),p3(2),p3(3),'ob')
plot3(te(1,:),te(2,:),te(3,:),'r','linewidth',2) 