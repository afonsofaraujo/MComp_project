function boundaryNodes = findNodesOnBoundary(fronteira, ymax, xval, coordx, coordy)
    indices = find(fronteira(:, 1) == xval);
    boundaryNodes = [indices, fronteira(indices, :)];
end