function u = kuQP_explicit(lambda, calG, calF, LFV, LGV, V)
    a = LFV + lambda * V;
    if a > 0
        u = -(LGV * calG) \ (LFV + LGV*calF + lambda*V);
    else
        u = -calG \ calF;
    end
end