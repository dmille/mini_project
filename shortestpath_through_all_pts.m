function ordered_pts = shortestpath_through_all_pts(points)
% shortestpath_through_all_pts Shortest path through all points
% Given a set of points, outputs a reordering of the points
% such that traveling along the weighted path, point i to j when i+1=j, is
% minimal.

    [MST,weights_mat] = EMST(points);
    sp_order = MST_to_shortestpath(MST,weights_mat);
    ordered_pts = points(sp_order,:);
    
end

