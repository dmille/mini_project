function line = line_param(A, B)
%LINE_PARAM: line parameters (L): x + m*y + b = 0 with L is the line
%through A and B
    if(A(1) == B(1) && A(2) ~= B(2))
        line = [1 0 -A(1)];
    elseif(A(2) == B(2) && A(1) ~= B(1))
        line = [0 1 -A(2)];
    else
        L = [A(2) 1; B(2) 1];
        b = [-A(1) -B(1)];
        sol = L^-1*b';
        line = [1 sol'];
    end
end

