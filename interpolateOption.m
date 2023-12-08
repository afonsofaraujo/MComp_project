


% [fineXY, fineU, fineP] = interpolateOption(coordx, coordy, u, xcentroid, ycentroid, pressure, fronteira, 1);
% % Scatter plot for fineU
% figure;
% scatter(fineXY(:, 1), fineXY(:, 2), 10, fineU, 'filled');
% colorbar;
% title('Fine U Distribution');
% xlabel('X-axis');
% ylabel('Y-axis');
% axis equal;
% 
% % Scatter plot for fineP
% figure;
% scatter(fineXY(:, 1), fineXY(:, 2), 10, fineP, 'filled');
% colorbar;
% title('Fine P Distribution');
% xlabel('X-axis');
% ylabel('Y-axis');
% axis equal;


% % Plot boundary
% B5 = [B1(:,1), 2*max(coordy)-B1(:,2)];
% B6 = [B3(:,1), 2*max(coordy)-B3(:,2)];
% B7 = [B2(:,1), 2*max(coordy)-B2(:,2)];
% plot(B1(:,1), B1(:,2),'k', 'LineWidth', 1),
% plot(B2(:,1), B2(:,2),'k', 'LineWidth', 1);
% plot(B3(:,1), B3(:,2),'k', 'LineWidth', 1);
% plot(B5(:,1), B5(:,2),'k', 'LineWidth', 1);
% plot(B6(:,1), B6(:,2),'k', 'LineWidth', 1);
% plot(B7(:,1), B7(:,2),'k', 'LineWidth', 1);




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
