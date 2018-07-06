function [out_mat,out_pts] = add_natural_ending(in_mat,interp_t,interp_p)
%add_natural_ending Natural ending to the b-sline in the coeff matrix before \
%   Detailed explanation goes here
    %defining natural end conditions for coeff matrix
    m = length(in_mat(1,:));
    
    s0 = interp_t(1);
    s1 = interp_t(2);
    s2 = interp_t(3);
    sn = interp_t(end);
    sn1 = interp_t(end-1);
    sn2 = interp_t(end-2);

    alpha0 = s2 - s0;
    gamma0 = s1 - s0;
    beta0 = -1*(alpha0) - gamma0;

    alpha_n = sn - sn1;
    gamma_n = sn - sn2;
    beta_n = -1*(alpha_n) - gamma_n;

    coeff_vec0 = zeros(1,m);
    coeff_vec_n = zeros(1,m);

    coeff_vec0(1:3) = [alpha0 beta0 gamma0];
    coeff_vec_n(end-2:end) = [alpha_n beta_n gamma_n];

    out_mat = [in_mat(1,:) ; coeff_vec0 ; in_mat(2:end-1,:) ; coeff_vec_n ; in_mat(end,:)];
    
    %defining natural ending conditions for interpolation points
    out_pts = [interp_p(1,:) ; 0 0 ; interp_p(2:end-1,:) ; 0 0 ; interp_p(end,:)];
end

