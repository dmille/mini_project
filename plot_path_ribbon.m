function [] = plot_path_ribbon(pv, curve, width, norm_offset)

up = [0,0,1];

vectors = circshift(curve,-1) - curve;

vectors = [vectors(1:end-1, :); vectors(end-1,:)];
norm_vectors = vectors ./ vecnorm(vectors,2,2);

m = size(norm_vectors,1);

ups = repmat(up, m, 1);

left = cross(ups, norm_vectors, 2);
%left = left ./ vecnorm(left,2,2);

normals = cross(norm_vectors, left, 2);

left = left * width/2;
right = -left;

X = zeros(m, 3);
Y = zeros(m, 3);
Z = zeros(m, 3);

X(:,2) = curve(:,1) + normals(:,1) * norm_offset;
Y(:,2) = curve(:,2) + normals(:,2) * norm_offset;
Z(:,2) = curve(:,3) + normals(:,3) * norm_offset;

X(:,1) = X(:,2) - left(:,1);
Y(:,1) = Y(:,2) - left(:,2);
Z(:,1) = Z(:,2) - left(:,3);

X(:,3) = X(:,2) - right(:,1);
Y(:,3) = Y(:,2) - right(:,2);
Z(:,3) = Z(:,2) - right(:,3);

surf(pv, X,Y,Z);
%plot3(pv, curve(:,1), curve(:,2), curve(:,3)-0.005, ':r');

end