function success = write_off(xyz, filename)

%% Prepare the data

    num_vertices = size(xyz,1);
    heightmap = reshape(xyz(:,3), [256,256]);
    [m n] = size(heightmap);

    xyz_scaled = rescale(xyz);
    num_faces = (m-1) * (n-1) * 2;


    %% Open File for writing
    fid = fopen(filename,'w');

    %% Write header with number of vertices and faces
    fprintf(fid, 'OFF\n');
    fprintf(fid, '%d %d 0\n', num_vertices, num_faces);

    %% Write out all the vertices
    for i = [1:num_vertices]       
        fprintf(fid,'%f %f %f\n', xyz_scaled(i,1), xyz_scaled(i,2), xyz_scaled(i,3));
    end

    index_of = @(i,j) (j-1)*m+i-1;

    %% Write out all the faces
    for i = [1:m-1]
        for j = [1:n-1]
            
            %% write lower triangle
            fprintf(fid, '3 %d %d %d\n', index_of(i,j), index_of(i+1, j), ...
                    index_of(i+1, j+1));
            
            %% write upper triangle
            fprintf(fid, '3 %d %d %d\n', index_of(i,j), index_of(i+1, j+1), ...
                    index_of(i, j+1));
            
        end
    end
        
    fclose(fid);
end