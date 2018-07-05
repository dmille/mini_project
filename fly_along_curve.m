function [] = fly_along_curve(curve, pv, mmv, speed)

    n_points = size(curve, 1);
    hlight = camlight(pv, 'headlight');
    delay = 1.0/speed;

    for i = [1:n_points-5]
        
        campos(pv, curve(i,:));
        camtarget(pv, curve(i+5,:));
        camlight(hlight,'headlight');
        drawnow;

        %set(pv, 'CameraTarget', test_curve(i+1,:));
        pause(delay);
    end



end