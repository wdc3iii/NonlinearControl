function v = kvQP(eta, lambda, P, F, G)
    v = quadprog(1, 0, 2*eta'*P*G, -lambda*eta'*P*eta - eta'*(F'*P + P*F)*eta);
end
