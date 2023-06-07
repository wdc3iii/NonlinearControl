function u = kQP(lambda, V, LfV, LgV)
    u = quadprog(1, 0, LgV, -LfV - lambda * V);
end
