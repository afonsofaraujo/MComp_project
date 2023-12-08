function [Ke, fe]=Elem_Quad4 (XN,fL)
    Ke = zeros(4,4);
    fe = zeros(4,1);
    nip = 9;
    [xp, wp]=Genip2DQ (nip);
    for ip=1:nip	%	ciclo para os pontos de integracao
        %seleção do ponto a tratar
        csi = xp(ip,1);
        eta = xp(ip,2);
        %transformação de coordenadas
        [B, psi, Detj]=Shape_N_Der4 (XN,csi,eta);
        %adição da contr
        wip = wp(ip)*Detj;
        Ke = Ke + wip*B*B' ;    %   matriz de rigidez
        wipf = fL*wip ;
        fe = fe + wipf*psi ;    %   vector de forcas (nulo no nosso caso)
        disp(fe)
    end
end