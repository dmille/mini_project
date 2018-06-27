for i=1:255
    for j=1:255
        faceU = ['3 ',num2str(((i-1)*256)+j),' ',num2str(((i-1)*256)+j+1),' ',num2str(((i-1)*256)+j+257)];
        disp(faceU);
        faceL = ['3 ',num2str(((i-1)*256)+j),' ',num2str(((i-1)*256)+j+256),' ',num2str(((i-1)*256)+j+257)];
        disp(faceL);
    end
end