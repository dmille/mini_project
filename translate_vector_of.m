function [transpose_vector] = translate_vector_of(A, B)
%TRANSPOSE_VECTOR_OF Summary of this function goes here
%   Detailed explanation goes here
    transpose_vector = [B(1) - A(1), B(2) - A(2)];
end

