function [sp_order] = MST_to_shortestpath(MST,weights_mat)
% MST_to_shortestpath Computes shortest path order on a set of nodes in a weighted minimum
% spanning tree. distmat is a matrix containing all the weights for a
% completely connected graph with the nodes of the MST.
%   Detailed explanation goes here

    while max(degree(MST)) > 2
        deg = degree(MST);
        deg3_ind = find(deg > 2);
        deg1_ind = find(deg == 1);

        num_deg3 = length(deg3_ind);
        num_deg1 = length(deg1_ind);

        add_del_mat = zeros(num_deg3*num_deg1*(num_deg1-1)*2,5);
        for k=1:num_deg3    
            for i=1:num_deg1
                sp = shortestpath(MST,deg3_ind(k),deg1_ind(i));
                for j = 1:num_deg1
                    if (i ~= j)
                        if (~ismember([min(sp(end),deg1_ind(j)) max(sp(end),deg1_ind(j))],MST.Edges.EndNodes,'rows'))
                            add_del_mat(find(add_del_mat(:,1) == 0,1),:) = [sp(1) sp(2) sp(end) deg1_ind(j) weights_mat(sp(end),deg1_ind(j))-weights_mat(sp(1),sp(2))];
                        end
                        if (~ismember([min(sp(2),deg1_ind(j)) max(sp(2),deg1_ind(j))],MST.Edges.EndNodes,'rows'))
                            add_del_mat(find(add_del_mat(:,1) == 0,1),:) = [sp(1) sp(2) sp(2) deg1_ind(j) weights_mat(sp(2),deg1_ind(j))-weights_mat(sp(1),sp(2))];
                        end
                    end
                end
            end
        end

        %min_add_pos_dist = min(add_del_mat(add_del_mat(:,5)>0,5));
        add_del_mat = add_del_mat(add_del_mat(:,1) ~= 0,:);
        add_del_order = sortrows(add_del_mat,5);
        %min_row = add_del_mat(find(add_del_mat(:,5) == min_add_pos_dist,1),:);
        inter_MST = MST;
        n=1;
        min_row = add_del_order(n,:);
        inter_MST = rmedge(inter_MST,min(min_row(1:2)),max(min_row(1:2)));
        inter_MST = addedge(inter_MST,min(min_row(3:4)),max(min_row(3:4)),weights_mat(min(min_row(3:4)),max(min_row(3:4))));

        while max(conncomp(inter_MST) > 1)
            inter_MST = MST;
            n = n+1;
            min_row = add_del_order(n,:);
            inter_MST = rmedge(inter_MST,min(min_row(1:2)),max(min_row(1:2)));
            inter_MST = addedge(inter_MST,min(min_row(3:4)),max(min_row(3:4)),weights_mat(min(min_row(3:4)),max(min_row(3:4))));
        end

        MST = inter_MST;

    end

    path_end_points = find(degree(MST) == 1);
    sp_order = shortestpath(MST,path_end_points(1),path_end_points(2));
    
end

