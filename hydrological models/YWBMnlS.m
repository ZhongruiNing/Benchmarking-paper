function [Q, E, SN] = YWBMnlS(x, params)

    P   = x(:, 1);
    T   = x(:, 2);
    PET = x(:, 3);
    
    Ks      = params(1);
    Kg      = params(2);
    Ksn     = params(3);
    alpha   = params(4);
    smax    = params(5);
    
    n   = length(P);
    Pr  = zeros(n, 1);    % Rainfall
    Psn = zeros(n, 1);    % Snowfall
    
    Qs  = zeros(n, 1);    % Surface runoff
    Qg  = zeros(n, 1);    % Underground runoff
    Q   = zeros(n, 1);    % Total runoff
    
    S   = zeros(n, 1);    % Soil moisture
    SN  = zeros(n, 1);    % Snowpack
    E   = zeros(n, 1);    % Actual evapotranspiration
    
    Tm  = -2;           % Lower threshold of snow
    Ts  = 4;            % Upper threshold of snow
    
    S(1)  = smax / 3;     % Initial soil moisture
    SN(1) = 0;            % Initial snowpack
    
    
    
    % Main body for WBM
    for i = 2 : n
        S_p  = S(i - 1);
        SN_p = SN(i - 1);
    
        % Snow module
        % if min(T) < 0
        snow_ratio = max([1 - exp(-1 * min((T(i) - Ts) / (Ts - Tm), 0)^2), 0]);

        % Separate precipitation to rainfall and snowfall
        Psn(i) = P(i) * snow_ratio;
        Pr(i) = P(i) * (1 - snow_ratio);

        % Snowfall accumulates to snowpack
        snow = SN_p + Psn(i);

        % Snow melting ratio
        melt_ratio = min([max([1 - exp(-min((Tm - T(i) - 4) / (Ts - Tm + 1), 0)^6), 0]), 1]);


        % Snow melting processing
        melt = Ksn * snow * melt_ratio;
        snow = snow - melt;
        SN(i) = snow;
        Pr(i) = Pr(i) + melt;

        % Underground runoff based on runoff in last time step
        Qg(i) = Kg * S_p;
        water_available = S_p - Qg(i) + Pr(i);
        
        % Evaporation Module
        temp_x = PET(i) / water_available;
        w = 1 / (1 - alpha);
        actual_evap = water_available * (1 + temp_x - (1 + temp_x ^ w) ^ (1 / w));
        E(i) = E(i) + actual_evap;

        water_available = water_available - actual_evap;

        water_available = max(0.001, water_available);

        % Runoff generation module
        % Surface runoff
        Qs1 = max(0, water_available - smax);
        water_available = water_available - Qs1;
        Qs2 = Ks * water_available * tanh(Pr(i) / smax);
        Qs(i) = Qs1 + Qs2;

        water_available = water_available - Qs2;

        S(i) = water_available;
        
        Q(i) = Qs(i) + Qg(i);
    end

end
