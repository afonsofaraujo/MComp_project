function [Kg, fg] = applyBoundaryConditions(Kg, fg, B1, B2, B4, coordout, U)
    boom = 1.0e+14;
    for i = 1:length(B1)
        j = find(coordout(:, 2) == B1(i, 1) & coordout(:, 3) == B1(i, 2));
        Kg(j, j) = boom;
        fg(j) = boom * coordout(j, 3) * U;
    end

    for i = 1:length(B2)
        j = find(coordout(:, 2) == B2(i, 1) & coordout(:, 3) == B2(i, 2));
        Kg(j, j) = boom;
        fg(j) = 0;
    end

    for i = 1:length(B4)
        j = find(coordout(:, 2) == B4(i, 1) & coordout(:, 3) == B4(i, 2));
        Kg(j, j) = boom;
        fg(j) = boom * max(B4(:, 2)) * U;
    end

    disp('Boundary conditions applied...');
end