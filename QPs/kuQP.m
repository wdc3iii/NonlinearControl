function u = kuQP(lambda, calG, calF, LFV, LGV, V)
    u = quadprog(calG'*calG, calF'*calG, LGV*calG, -lambda*V - LFV - LGV*calF);
end