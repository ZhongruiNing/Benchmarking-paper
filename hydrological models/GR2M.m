function [Q, S2] = GR2M(x, params)

    x1 = params(1);
    x2 = params(2);

    P = x(:, 1);
    E = x(:, 3);

    n = length(P);

    S = nan(n, 1);
    Q = nan(n, 1);
    R = nan(n, 1);
    S2 = nan(n, 1);

    S0 = x1;
    R0 = 0;

    for i = 1 : n
        vphi = tanh(P(i) / x1);
        psi  = tanh(E(i) / x1);

        S1   = (S0 + x1 * vphi) / (1 + vphi * (S0 / x1));
        P1   = P(i) + S0 - S1;
        S2(i)   = S1 * (1 - psi) / (1 + psi * (1 - (S1 / x1)));

        S(i) = S2(i) / ((1 + (S2(i) / x1) ^ 3) ^ (1/3));

        P2   = S2(i) - S0;
        P3   = P1 + P2;
        R1   = R0 + P3;
        R2   = x2 * R1;

        Q(i) = (R2 ^ 2) / (R2 + 60);
        R(i) = R2 - Q(i);

        S0   = S(i);
        R0   = R(i);
    end

end
