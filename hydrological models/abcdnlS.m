function [Q, Ea, SN] = abcdnlS(x, params)

    a   = params(1);
    b   = params(2);
    c   = params(3);
    d   = params(4);
    Ksn = params(5);

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

    Tm  = -2;           % Lower threshold of snow
    Ts  = 4;            % Upper threshold of snow

    for i = 2 : n
        SN_P = SN(i - 1);
        % Snow module
        % if min(T) < 0
        snow_ratio = max([1 - exp(-1 * min((T(i) - Ts) / (Ts - Tm), 0)^2), 0]);

        % Separate precipitation to rainfall and snowfall
        Psn = P(i) * snow_ratio;
        Pr = P(i) * (1 - snow_ratio);

        % Snowfall accumulates to snowpack
        snow = SN_P + Psn;

        % Snow melting ratio
        melt_ratio = min([max([1 - exp(-min((Tm - T(i) - 4) / (Ts - Tm + 1), 0)^6), 0]), 1]);

        % Snow melting processing
        melt = Ksn * snow * melt_ratio;
        snow = snow - melt;
        SN(i) = snow;
        Pr = Pr + melt;

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
