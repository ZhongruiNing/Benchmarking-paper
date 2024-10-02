function [Q, AET, SN] = DWBMlS(x, params)
    alpha1  = params(1);
    alpha2  = params(2);
    smax    = params(3);
    d       = params(4);
    m       = params(5);

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
    AET     = zeros(n, 1);
    S       = zeros(n, 1);
    G       = zeros(n, 1);
    Pr      = zeros(n, 1);
    Ps      = zeros(n, 1);
    SN      = zeros(n, 1);

    Train = 4;
    Tsnow = -4;

    P(P == 0) = 0.01;
    for i = 2 : n
        SN_P = SN(i - 1);
        if T(i) < Tsnow
            Ps(i) = P(i);
            Pr(i) = 0;
        elseif T(i) > Train
            Ps(i) = 0;
            Pr(i) = P(i);
        else
            Ps(i) = P(i) * (Train - T(i)) / (Train - Tsnow);
            Pr(i) = P(i) - Ps(i);
        end
        SN(i) = SN_P + Ps(i);

        if T(i) < Tsnow
            SM = 0;
        elseif T(i) > Train
            SM = m * SN(i);
        else
            SM = m * SN(i) * (Train - T(i)) / (Train - Tsnow);
        end
        SN(i) = SN(i) - SM;
        Pr(i) = Pr(i) + SM;
        
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
        check_complex(Y(i))

        R(i)   = W(i) - Y(i);
        AET(i) = W(i) * Fu(PET(i) / W(i), alpha2);

        S(i)   = Y(i) - AET(i);

        Qb(i)  = d * G0;
        G(i)   = (1 - d) * G0 + R(i);

        Q(i)   = Qd(i) + Qb(i);
        check_complex(Q(i))
        S0     = S(i);
        check_complex(S0);
        G0     = G(i);
    end
end