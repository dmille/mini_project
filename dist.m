function distance = dist(A,B)
%DIST: calculate distance between 2 points (2-D)
    distance = pdist([A(1), A(2); B(1),B(2)], 'euclidean');
end

