function [] = fly_along_curve(curve, pv, mmv, speed, minimap_cam_height)

    n_points = size(curve, 1);
    hlight = camlight(pv, 'headlight');
    delay = 1.0/speed;

    plane_marker = plot3(mmv,curve(1,1),curve(1,2),curve(1,3),'b.', ...
                         'MarkerSize', 20);

    for i = [1:n_points-5]
        set(plane_marker, 'Xdata', curve(i,1));
        set(plane_marker, 'Ydata', curve(i,2));
        set(plane_marker, 'Zdata', curve(i,3));

        campos(mmv, [curve(i,1:2), minimap_cam_height]);
        camtarget(mmv, [curve(i,1:2), 0]);

        campos(pv, curve(i,:));               
        camtarget(pv, curve(i+5,:));

        camlight(hlight,'headlight');
        drawnow;

        pause(delay);
    end



end