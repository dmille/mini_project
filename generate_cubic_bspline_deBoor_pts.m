function [d_pts] = generate_cubic_bspline_deBoor_pts(interp_p)
%generate_cubic_bspline_deBoor_pts given a set of points, interp_p, output
%is a set of deBoor points for a B-spline that passes through all interp_p
% on the time interval 0-1. 
%   Detailed explanation goes here
    
    % cubic B-spline
    k = 4; 
    
    % gets cord length t for generating the coeff matrix and the t knots
    % sequence (duplicate first and last)
    [interp_t,t_knot] = chord_length_parametrization(interp_p); 
      
    % gets the coeff matrix to be used in M\p
    coeff_mat = generate_cubic_bspline_coefficent_matrix(t_knot, interp_t);
    
    % addes natural ending in order to solve system of equations
    [coeff_mat,interp_p] = add_natural_ending(coeff_mat,interp_t,interp_p);
    
    % solves system for deBoor points
    d_pts = coeff_mat\interp_p;
end

