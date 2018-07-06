function [out_mat] = generate_cubic_bspline_coefficent_matrix(t_knot,interp_t)
%UNTITLED6 Given the tknots, generates a coefficent matrix to be used to as
% M*D = P, Where M is the output matrix, D are the deBoor control points of the
% spline, and P are the points on the spline. 
%   Detailed explanation goes here
    k = 4;
    n = length(t_knot) - k - 1;
    out_mat = zeros(length(interp_t),n+1);
    for i=1:length(interp_t)
        t = interp_t(i);

        %t = 4;
        if (t == t_knot(end))
            r = length(t_knot) - k;
        else
            r = find(t_knot > t,1) - 1;
        end

        a20 = (t_knot(r+1) - t)/(t_knot(r+1) - t_knot(r));
        a11 = (t_knot(r+2) - t)/(t_knot(r+2) - t_knot(r));
        a10 = (t_knot(r+1) - t)/(t_knot(r+1) - t_knot(r-1));
        a02 = (t_knot(r+3) - t)/(t_knot(r+3) - t_knot(r));
        a01 = (t_knot(r+2) - t)/(t_knot(r+2) - t_knot(r-1));
        a00 = (t_knot(r+1) - t)/(t_knot(r+1) - t_knot(r-2));

        c1 = a20 * a10 * a00;
        c2 = a20 * a10 * (1 - a00) + (a20 * (1 - a10) + a11 * (1 - a20)) * a01;
        c3 = (1 - a20) * (1 - a11) * a02 + (a20 * (1 - a10) + a11 * (1 - a20)) * (1 - a01);
        c4 = (1 - a20) * (1 - a11) * (1 - a02);
        coeff = [c1 c2 c3 c4];

        coeff_vec = zeros(1,n+1);
        coeff_vec(r-k+1:r) = coeff;

        out_mat(i,:) = coeff_vec;
    end

end

