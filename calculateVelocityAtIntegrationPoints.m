function [Res, xint, yint, u1, v1] = calculateVelocityAtIntegrationPoints(connectivityData, coordx, coordy, u, nip, elementType)
    Nels = size(connectivityData, 1);
    Res = 0;
    xint = [];
    yint = [];
    arrow_u = [];
    arrow_v = [];
    rotationMatrix = [cosd(90) -sind(90); sind(90) cosd(90)];
    for i = 1:Nels
        edofs = connectivityData(i,:);
        XN = [coordx(edofs), coordy(edofs)];
        [xp, wp] = Genip2DQ(nip);
        for ip = 1:nip
            csi = xp(ip, 1);
            eta = xp(ip, 2);
            if strcmpi(elementType, 'QUAD4')
                [B, psi, Detj] = Shape_N_Der4(XN, csi, eta);
            elseif strcmpi(elementType, 'QUAD8')
                [B, psi, Detj] = Shape_N_Der8(XN, csi, eta);
            else
                error('Invalid element type. Supported types are Q4 and Q8.');
            end
            uint = psi' * u(edofs);
            Res = Res + uint * Detj * wp(ip);
            xpint = XN' * psi;
            gradu = B' * u(edofs);
            fluxu = -gradu;
            xint = [xint; xpint(1)];
            yint = [yint; xpint(2)];
            arrow_u = [arrow_u; fluxu(1)];
            arrow_v = [arrow_v; fluxu(2)];
        end
    end
    rotatedVectors = rotationMatrix * [arrow_u(:)'; arrow_v(:)'];
    u1 = reshape(rotatedVectors(1, :), size(xint));
    v1 = reshape(rotatedVectors(2, :), size(yint));
end