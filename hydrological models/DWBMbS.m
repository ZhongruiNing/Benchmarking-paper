function [Q, AET, SN] = DWBMbS(x, params)
    alpha1  = params(1);
    alpha2  = params(2);
    smax    = params(3);
    d       = params(4);
    beta    = params(5);

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
        SN0 = SN(i - 1);
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
        RFt = Pr(i);
        SFt = Ps(i);
        if T(i) < Tsnow
            SN0 = SN0 + S0 - 0.1;
            S0  = 0.1;
        end
        a = (RFt + S0) / (RFt + SFt + S0 + SN0);
        b = (SFt + SN0) / (RFt + SFt + S0 + SN0);
        
        
        PET(i) = PET(i) * a;
        if PET(i) == 0
            PET(i) = 0.1;
        end
        PMt = PET(i) * (b/0.15);

        if PMt == 0
            Mt = 0;
        else
            Mt  = (SFt + SN0) * Fu(PMt / (SFt + SN0), beta);
        end
        if isnan(Mt)
            Mt = 0;
        end
        check_complex(Mt);

        if isnan(Mt)
            Mt = 0;
        end

        SNt = SN0 + SFt - Mt;
        if SNt < 0
            SNt = 0;
        end
        
        SN(i)  = SNt;
        P(i)   = RFt + Mt;

        if T(i) < Tsnow
            P(i) = max(0.001, P(i));
            X(i) = 0;
            Qd(i) = P(i) - X(i);
        else
            X0     = smax - S0 + PET(i);
            P(i) = max(0.001, P(i));
            X(i)   = P(i) * Fu(X0 / P(i), alpha1);
            if isinf(X(i))
                X(i) = P(i) - 0.0001;
            end

            Qd(i)  = P(i) - X(i);
            Qd(i)  = max([0, Qd(i)]);
        end
        
        W(i)   = X(i) + S0;
        Y0     = smax + PET(i);

        Y(i)   = W(i) * Fu(Y0 / W(i), alpha2);
        check_complex(Y(i));

        R(i)   = W(i) - Y(i);
        AET(i) = W(i) * Fu(PET(i) / W(i), alpha2);
        check_complex(AET(i));

        S(i)   = Y(i) - AET(i);


        Qb(i)  = d * G0;
        G(i)   = (1 - d) * G0 + R(i);

        Q(i)   = Qd(i) + Qb(i);
        S0     = S(i);
        S0     = max(0.1, S0);
        G0     = G(i);

        
    end
end