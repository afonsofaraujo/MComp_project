function [fineXY, fineU, fineP] = interpolateOption(coordx, coordy, u, xcentroid, ycentroid, pressure, fronteira, fineMeshDensity)
    fineX = linspace(min(coordx), max(coordx), fineMeshDensity * (length(coordx) - 1) + 1);
    fineY = linspace(min(coordy), max(coordy), fineMeshDensity * (length(coordy) - 1) + 1);
    [fineX, fineY] = meshgrid(fineX, fineY);
    fineU = griddata(coordx, coordy, u, fineX, fineY);
    fineP = griddata(xcentroid, ycentroid, pressure, fineX, fineY);
    inBoundary = inpolygon(fineX, fineY, fronteira(:, 1), fronteira(:, 2));
    fineU(~inBoundary) = NaN;
    fineP(~inBoundary) = NaN;

    % Reshape into single columns and combine coordinates
    fineXY = [fineX(:), fineY(:)];
    fineU = fineU(:);
    fineP = fineP(:);
end
