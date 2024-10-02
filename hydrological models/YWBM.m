function [Q, E] = YWBM(x, params)

    P   = x(:, 1);
    PET = x(:, 3);
    
    Ks      = params(1);
    Kg      = params(2);
    alpha   = params(3);
    smax    = params(4);
    
    n   = length(P);
    
    Qs  = zeros(n, 1);    % Surface runoff
    Qg  = zeros(n, 1);    % Underground runoff
    Q   = zeros(n, 1);    % Total runoff
    
    S   = zeros(n, 1);    % Soil moisture
    E   = zeros(n, 1);    % Actual evapotranspiration
    
    S(1)  = smax / 3;     % Initial soil moisture
    
    % Main body for WBM
    for i = 2 : n
        S_p  = S(i - 1);

        % Underground runoff based on runoff in last time step
        Qg(i) = Kg * S_p;
        water_available = S_p - Qg(i) + P(i);
        
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
        Qs2 = Ks * water_available * tanh(P(i) / smax);
        Qs(i) = Qs1 + Qs2;

        water_available = water_available - Qs2;

        S(i) = water_available;
        
        Q(i) = Qs(i) + Qg(i);
    end
end
