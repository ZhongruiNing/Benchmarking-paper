function [Q, Ea, SN] = abcdbS(x, params)

    a   = params(1);
    b   = params(2);
    c   = params(3);
    d   = params(4);
    beta= params(5);

    P   = x(:, 1);
    T   = x(:, 2);
    PET = x(:, 3);
    Inv = [200, 5];

    S0  = Inv(1);
    S   = Inv(2);
    SWE = S0 + S;

    n   = length(P);

    Q   = zeros(n, 1);
    Q1  = zeros(n, 1);
    Q2  = zeros(n, 1);
    SN  = zeros(n, 1);
    Ea  = zeros(n, 1);

    Pr      = zeros(n, 1);
    Ps      = zeros(n, 1);

    Train = 0.5;
    Tsnow = -0.5;

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
        aa = (RFt + S0) / (RFt + SFt + S0 + SN0);
        bb = (SFt + SN0) / (RFt + SFt + S0 + SN0);
      
        
        PET(i) = PET(i) * aa;
        if PET(i) == 0
            PET(i) = 0.1;
        end
        PMt = PET(i) * (bb/0.15);

        if PMt == 0
            Mt = 0;
        else
            Mt  = (SFt + SN0) * Fu(PMt / (SFt + SN0), beta);
        end
        if isnan(Mt)
            Mt = 0;
        end
        check_complex(Mt);
        
        SN(i) = SN0 + SFt - Mt;
        
        P(i) = RFt + Mt;

        if T(i) > Tsnow
            W = P(i) + S0;
            Y = 0.5 * (W + b) / a - ((0.5 * (W + b) / a) .^ 2 - W * b / a) .^ 0.5;
            S0 = Y * exp(-PET(i) / b);

            Ea(i) = min(PET(i), max(0, Y * (1 - exp(-PET(i) / b))));
            
            %S=S+c*(W-Ea(i));
            S = (c * (W - Y) + S) / (1 + d);
            Q1(i) = (1 - c) * (W - Y);
            Q2(i) = d * S;
            Q(i) =  Q1(i) + Q2(i);
            S = (1 - d) * S;
            SWE(i) = S0 + S;
        else
            Q1(i) = max(P(i) - 0.01, 0);
            P(i) = min(P(i), 0.01);
            W = P(i) + S0;
            Y = 0.5 * (W + b) / a - ((0.5 * (W + b) / a) .^ 2 - W * b / a) .^ 0.5;
            S0 = Y * exp(-PET(i) / b);

            Ea(i) = min(PET(i), max(0, Y * (1 - exp(-PET(i) / b))));

            S = (c * (W - Y) + S) / (1 + d);
            Q1(i) = Q1(i) + (1 - c) * (W - Y);
            Q2(i) = d * S;
            Q(i) =  Q1(i) + Q2(i);
            S = (1 - d) * S;
            SWE(i) = S0 + S;
        end
    end

end
