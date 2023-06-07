function u = ktQP_explicit(lambda, yd, calG, calF, LFV, LGV, V)
    a = LFV + lambda * V;
    if a > 0
        u = -(LGV * calG) \ (LFV + LGV*(calF-yd) + lambda*V);
    else
        u = -calG \ (calF - yd);
    end
end