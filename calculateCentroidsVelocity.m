% function [vx, vy] = calculateCentroidsVelocity(connectivityData, coordx, coordy, u)
%     Nels = size(connectivityData, 1);
%     fluxArray = zeros(Nels, 2);
%     rotationMatrix = [cosd(90) -sind(90); sind(90) cosd(90)];
%     for i = 1:Nels
%         edofs = connectivityData(i, :);
%         XN = [coordx(edofs), coordy(edofs)];
%         [B, ~, ~] = Shape_N_Der4(XN, 0, 0);
%         gradu = B' * u(edofs);
%         fluxu = -gradu;
%         fluxArray(i, :) = fluxu;
%     end
%     rotatedVectors = rotationMatrix * [fluxArray(:, 1)'; fluxArray(:, 2)'];     % Rotate flux vectors
%     vx = rotatedVectors(1, :)';
%     vy = rotatedVectors(2, :)';
% end
function [vx, vy] = calculateCentroidsVelocity(connectivityData, coordx, coordy, u, elementType)
    Nels = size(connectivityData, 1);
    fluxArray = zeros(Nels, 2);
    rotationMatrix = [cosd(90) -sind(90); sind(90) cosd(90)];
    for i = 1:Nels
        edofs = connectivityData(i, :);
        XN = [coordx(edofs), coordy(edofs)];

        if strcmpi(elementType, 'QUAD4')
            [B, ~, ~] = Shape_N_Der4(XN, 0, 0);
        elseif strcmpi(elementType, 'QUAD8')
            [B, ~, ~] = Shape_N_Der8(XN, 0, 0); % Replace with the actual function for Q8
        else
            error('Invalid element type. Supported types are Q4 and Q8.');
        end

        gradu = B' * u(edofs);
        fluxu = -gradu;
        fluxArray(i, :) = fluxu;
    end
    rotatedVectors = rotationMatrix * [fluxArray(:, 1)'; fluxArray(:, 2)'];     % Rotate flux vectors
    vx = rotatedVectors(1, :)';
    vy = rotatedVectors(2, :)';
end
