function [Q, AET] = DWBM(x, params)
    alpha1  = params(1);
    alpha2  = params(2);
    smax    = params(3);
    d       = params(4);

    P       = x(:, 1);
    PET     = x(:, 3);

    S0      = 0.5 * smax;
    G0      = 0.5 * smax;

    n       = length(P);
    Q       = nan(n, 1);
    Qd      = nan(n, 1);
    Qb      = nan(n, 1);
    X       = nan(n, 1);
    Y       = nan(n, 1);
    R       = nan(n, 1);
    W       = nan(n, 1);
    AET     = nan(n, 1);
    S       = nan(n, 1);
    G       = nan(n, 1);

    for i = 1 : n
        X0     = smax - S0 + PET(i);
        P(i) = max(0.001, P(i));
        X(i)   = P(i) * Fu(X0 / P(i), alpha1);
        if isinf(X(i))
            X(i) = P(i) - 0.0001;
        end

        Qd(i)  = P(i) - X(i);
        Qd(i)  = max([0, Qd(i)]);

        W(i)   = X(i) + S0;
        Y0     = smax + PET(i);

        Y(i)   = W(i) * Fu(Y0 / W(i), alpha2);

        R(i)   = W(i) - Y(i);
        AET(i) = W(i) * Fu(PET(i) / W(i), alpha2);

        S(i)   = Y(i) - AET(i);

        Qb(i)  = d * G0;
        G(i)   = (1 - d) * G0 + R(i);

        Q(i)   = Qd(i) + Qb(i);
        S0     = S(i);
        G0     = G(i);
    end
end