function [u, xcentroid, ycentroid, vx, vy, xint, yint, vxint, vyint, pressure] = readOutput(outputfileName)
    fid = fopen(outputfileName, 'r');
    % Skip 3 lines
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    Nnds = fscanf(fid, '%d', 1);
    fgetl(fid);
    u = fscanf(fid, '%f', [1, Nnds])';
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    Nels = fscanf(fid, '%d', 1);
    fgetl(fid);
    aux1 = fscanf(fid, '%f, %f, %f, %f', [4, Nels])';
    xcentroid = aux1(:,1);
    ycentroid = aux1(:,2);
    vx = aux1(:,3);
    vy = aux1(:,4);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    Nint = fscanf(fid, '%d', 1);
    fgetl(fid);
    aux2 = fscanf(fid, '%f, %f, %f, %f', [4, Nint])';
    xint = aux2(:,1);
    yint = aux2(:,2);
    vxint = aux2(:,3);
    vyint = aux2(:,4);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    Np = fscanf(fid, '%d', 1);
    fgetl(fid);
    aux3 = fscanf(fid, '%f, %f, %f', [3, Np])';
    pressure = aux3(:,3);
    fclose(fid);
end