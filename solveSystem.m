function u = solveSystem(Kg, fg)
    u = Kg \ fg;
    disp('System solved...');
end