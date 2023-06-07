function v = kvQP_explicit(eta, lambda, P, F, G)
    a = eta'*(F'*P + P*F)*eta + lambda*eta'*P*eta;
    b = G'*P*eta;
    if a > 0
        v = -0.5 * a /(b'*b)*b;
    else
        v = 0;
    end
end