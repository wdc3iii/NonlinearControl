 function u = kQP_explicit(lambda, V, LfV, LgV)
    a = LfV + lambda * V;
    b = LgV;
    if a <= 0
        u = 0;
    else
        u = -a / (b*b') * b;
    end
end

