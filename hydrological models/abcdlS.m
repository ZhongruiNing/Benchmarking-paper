function [Q, Ea, SN] = abcdlS(x, params)

    a   = params(1);
    b   = params(2);
    c   = params(3);
    d   = params(4);
    m   = params(5);

    P   = x(:, 1);
    T   = x(:, 2);
    PET = x(:, 3);
    Inv = [200, 5];

    H   = Inv(1);
    S   = Inv(2);
    SWE = H + S;

    n   = length(P);

    Q   = zeros(n, 1);
    Q1  = zeros(n, 1);
    Q2  = zeros(n, 1);
    SN  = zeros(n, 1);

    Train = 4;
    Tsnow = -4;
    


    for i = 2 : n
        SN_P = SN(i - 1);
        if T(i) < Tsnow
            Ps = P(i);
            Pr = 0;
        elseif T(i) > Train
            Ps = 0;
            Pr = P(i);
        else
            Ps = P(i) * (Train - T(i)) / (Train - Tsnow);
            Pr = P(i) - Ps;
        end
        SN(i) = SN_P + Ps;
        if T(i) < Tsnow
            SM = 0;
        elseif T(i) > Train
            SM = m * SN(i);
        else
            SM = m * SN(i) * (Train - T(i)) / (Train - Tsnow);
        end
        SN(i) = SN(i) - SM;
        Pr = Pr + SM;

        W = Pr + H;
        Y = 0.5 * (W + b) / a - ((0.5 * (W + b) / a) .^ 2 - W * b / a) .^ 0.5;
        H = Y * exp(-PET(i) / b);

        Ea(i) = min(PET(i), max(0, Y * (1 - exp(-PET(i) / b))));
        
        %S=S+c*(W-Ea(i));
        S = (c * (W - Y) + S) / (1 + d);
        Q1(i) = (1 - c) * (W - Y);
        Q2(i) = d * S;
        Q(i) =  Q1(i) + Q2(i);
        S = (1 - d) * S;
        SWE(i) = H + S;
    end

end
