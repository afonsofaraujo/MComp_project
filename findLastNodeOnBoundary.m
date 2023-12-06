function lastNode = findLastNodeOnBoundary(B3)
    lastNode = [max(B3(:, 1)), min(B3(:, 2))];
end