function line = ortho_line_param(A, B)
%ORTHO_LINE_PARAM: parameters of line (L'): -m*x + y + b2 = 0 
%which is orthogonal to the line      (L ): x + m*y + b1 = 0 at A with (L) is the line through A and B

% type -- 1: horizontal -- 2: vetical -- 3 : normal line

    temp = line_param(A,B);
    
    if(A(1) == B(1) && A(2) ~= B(2))
        line = [0 1 -A(2) 1];
    elseif(A(2) == B(2) && A(1) ~= B(1))
        line = [1 0 -A(1) 2];
    else
        line = [-temp(2), temp(1), temp(2)*A(1) - A(2), 3];
    end
end

