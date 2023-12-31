function [Ke, fe]=Elem_Quad4 (XN,fL)
    Ke = zeros(4,4);
    fe = zeros(4,1);
    nip = 9;
    [xp, wp]=Genip2DQ (nip);
    for ip=1:nip
        csi = xp(ip,1);
        eta = xp(ip,2);
        [B, psi, Detj]=Shape_N_Der4 (XN,csi,eta);
        wip = wp(ip)*Detj;
        Ke = Ke + wip*B*B' ;
        wipf = fL*wip;
        fe = fe + wipf*psi ;
    end
end