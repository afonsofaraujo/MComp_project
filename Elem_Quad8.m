function [Ke, fe] = Elem_Quad8(XN, fL)
    Ke = zeros(8, 8);
    fe = zeros(8, 1);
    nip = 9;
    [xp, wp] = Genip2DQ(nip);
    for ip = 1:nip
        csi = xp(ip, 1);
        eta = xp(ip, 2);
        [B, psi, Detj] = Shape_N_Der8(XN, csi, eta);
        wip = wp(ip) * Detj;
        Ke = Ke + wip * B * B';  %   stiffness matrix
        wipf = fL * wip;
        fe = fe + wipf * psi;    %   force vector
    end
end