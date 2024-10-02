function [Q, E] = TPWB(x, params)
    c   = params(1);
    SC  = params(2);

    P   = x(:, 1);
    PET = x(:, 3);

    n = length(P);

    Q = nan(n, 1);
    E = nan(n, 1);

    S0 = SC;

    for i = 1 : n
        E(i) = c * P(i) * tanh(P(i) / PET(i));
        S0   = S0 + P(i) - E(i);
        Q(i) = S0 * tanh(S0 / SC);
        S0   = S0 - Q(i);
    end
end