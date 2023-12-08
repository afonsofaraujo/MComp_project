function [xcentroid, ycentroid] = computeCentroids(connectivityData, coordx, coordy, elementType)
    Nels = size(connectivityData, 1);
    xcentroid = zeros(Nels, 1);
    ycentroid = zeros(Nels, 1);
    for i = 1:Nels
        edofs = connectivityData(i, :);
        XN = [coordx(edofs), coordy(edofs)];
        % Calculate centroid coordinates
        if strcmpi(elementType, 'QUAD4')
            xcentroid(i) = mean(XN(:, 1));
            ycentroid(i) = mean(XN(:, 2));
        elseif strcmpi(elementType, 'QUAD8')
            % For QUAD8, consider the mid-side nodes as well
            xcentroid(i) = mean([XN(:, 1); XN(5:8, 1)]);
            ycentroid(i) = mean([XN(:, 2); XN(5:8, 2)]);
        else
            error('Invalid element type. Supported types are QUAD4 and QUAD8.');
        end
    end
end