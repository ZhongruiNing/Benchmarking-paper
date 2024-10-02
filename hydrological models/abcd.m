function [Q, Ea] = abcd(x, params)
    a   = params(1);
    b   = params(2);
    c   = params(3);
    d   = params(4);

    P   = x(:, 1);
    PET = x(:, 3);
    Inv = [200, 5];

    H   = Inv(1);
    S   = Inv(2);
    SWE = H + S;

    n   = length(P);

    Q   = nan(n, 1);

    for i = 2 : n
        W = P(i) + H;
        Y = 0.5 * (W + b) / a - ((0.5 * (W + b) / a) .^ 2 - W * b / a) .^ 0.5;
        H = Y * exp(-PET(i) / b);

        Ea(i) = min(PET(i), max(0, Y * (1 - exp(-PET(i) / b))));
        
        %S=S+c*(W-Ea(i));
        S = (c * (W - Y) + S) / (1 + d);
        Q(i) = (1 - c) * (W - Y) + d * S;
        S = (1 - d) * S;
        SWE(i) = H + S;
    end
end