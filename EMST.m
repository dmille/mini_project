function [MST,weights_mat] = EMST(pts)
%EMST computes the euclidean minimum spanning tree for a set of points
%   Output is an weighted undirected graph MST
    
    weights_mat = zeros(length(pts));
    for i=1:length(pts)
       diffvec = pts(i,:) - pts;
       distvec = vecnorm(diffvec');
       weights_mat(i,:) = distvec;
    end

    G = graph(weights_mat,'upper','OmitSelfLoops');
    MST = minspantree(G);

end

