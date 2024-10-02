function [Q, AET, SN] = DWBMnlS(x, params)
    alpha1  = params(1);
    alpha2  = params(2);
    smax    = params(3);
    d       = params(4);
    Ksn     = params(5);

    P       = x(:, 1);
    T       = x(:, 2);
    PET     = x(:, 3);

    S0      = 0.5 * smax;
    G0      = 0.5 * smax;

    n       = length(P);
    Q       = zeros(n, 1);
    Qd      = zeros(n, 1);
    Qb      = zeros(n, 1);
    X       = zeros(n, 1);
    Y       = zeros(n, 1);
    R       = zeros(n, 1);
    W       = zeros(n, 1);
    E       = zeros(n, 1);
    AET     = zeros(n, 1);
    S       = zeros(n, 1);
    G       = zeros(n, 1);
    Pr      = zeros(n, 1);
    Ps      = zeros(n, 1);
    SN      = zeros(n, 1);

    Tm  = -2;           % Lower threshold of snow
    Ts  = 4;            % Upper threshold of snow

    P(P == 0) = 0.01;
    for i = 2 : n
        SN_P = SN(i - 1);
        % if min(T) < 0
        snow_ratio = max([1 - exp(-1 * min((T(i) - Ts) / (Ts - Tm), 0)^2), 0]);

        % Separate precipitation to rainfall and snowfall
        Ps(i) = P(i) * snow_ratio;
        Pr(i) = P(i) * (1 - snow_ratio);

        % Snowfall accumulates to snowpack
        SN(i) = SN_P + Ps(i);

        % SN(i) melting ratio
        melt_ratio = min([max([1 - exp(-min((Tm - T(i) - 4) / (Ts - Tm + 1), 0)^6), 0]), 1]);

        % SN(i) melting processing
        melt = Ksn * SN(i) * melt_ratio;
        SN(i) = SN(i) - melt;

        Pr(i) = Pr(i) + melt;
        
        P(i) = Pr(i);
        P(i) = max(0.001, P(i));

        X0     = smax - S0 + PET(i);
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