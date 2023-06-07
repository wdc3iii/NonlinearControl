function u = kSon(lambda, V, LfV, LgV)
    a = LfV + lambda * V;
    b = LgV';
    if b*b' == 0 
        u = 0;
    else
        u = -(a + sqrt(a^2 + (b*b')^2)) / (b*b') * b;
    end
end

