% interp_p are the points we are trying to interpolate using a B-spline. In
% this case we pass through each point in order.

interp_p = [-8 -8;-7 -7;-4 -6;-6 -4;-3 0;-2 4;1 7;2 9;5 0;6 -4;2 -8];

% d are the deBoor points we calculate based off the points we need to pass
% through
d = generate_cubic_bspline_deBoor_pts(interp_p);

% gets the tknots using chord length and tacking on leading and trailing
% duplicates. 
[interp_t,t_knot] = chord_length_parametrization(interp_p);

% sample t for plotting
tsampling = linspace(0,1,100);

% get bspline using tsampling
m = generate_cubic_bspline_coefficent_matrix(t_knot,tsampling);
bspline = m*d;

%plot
hold on;
plot(interp_p(:,1),interp_p(:,2),'*');
plot(d(:,1),d(:,2),'*');
plot(bspline(:,1),bspline(:,2));

