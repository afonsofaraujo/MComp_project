function pressure = calculatePressure(connectivityData, coordx, coordy, u, elementType, materialProperty)
    Nels = size(connectivityData, 1);
    constant = 104450;  % 101325+0.5*2.5^2*1000
    pressure = zeros(Nels, 1);
    fluxArray = zeros(Nels, 2); % Preallocate arrays for storing calculated values
    for i = 1:Nels
        edofs = connectivityData(i, :);
        XN = [coordx(edofs), coordy(edofs)];
        csi = 0;
        eta = 0;
        if strcmpi(elementType, 'QUAD4')
            [B, psi, Detj] = Shape_N_Der4(XN, csi, eta);
        elseif strcmpi(elementType, 'QUAD8')
            [B, psi, Detj] = Shape_N_Der8(XN, csi, eta);
        else
            error('Invalid element type. Supported types are QUAD4 and QUAD8.');
        end
        uint = psi' * u(edofs);
        gradu = B' * u(edofs);
        fluxu = -gradu;
        fluxArray(i, :) = fluxu;
        pressure(i) = (constant - 0.5*materialProperty*(sqrt(fluxu(1)^2 + fluxu(2)^2))^2) * 10^(-5);  % pressure in bar
    end
end