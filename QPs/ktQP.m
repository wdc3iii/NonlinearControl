function u = ktQP(lambda, yd, calG, calF, LFV, LGV, V)
    u = quadprog(calG'*calG, (calF-yd)'*calG, LGV*calG, -lambda*V - LFV - LGV*(calF - yd));
end