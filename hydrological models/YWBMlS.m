function [Q, E, SN] = YWBMlS(x, params)

    P   = x(:, 1);
    T   = x(:, 2);
    PET = x(:, 3);
    
    Ks      = params(1);
    Kg      = params(2);
    m       = params(3);
    alpha   = params(4);
    smax    = params(5);
    
    n   = length(P);
    
    Qs  = zeros(n, 1);    % Surface runoff
    Qg  = zeros(n, 1);    % Underground runoff
    Q   = zeros(n, 1);    % Total runoff
    
    S   = zeros(n, 1);    % Soil moisture
    SN  = zeros(n, 1);    % Snowpack
    E   = zeros(n, 1);    % Actual evapotranspiration
    
    Train = 4;
    Tsnow = -4;            % Upper threshold of snow
    
    S(1)  = smax / 3;     % Initial soil moisture
    SN(1) = 0;            % Initial snowpack
    
    
    
    % Main body for WBM
    for i = 2 : n
        S_p  = S(i - 1);
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

        % Underground runoff based on runoff in last time step
        Qg(i) = Kg * S_p;
        water_available = S_p - Qg(i) + Pr;
        
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
        Qs2 = Ks * water_available * tanh(Pr / smax);
        Qs(i) = Qs1 + Qs2;

        water_available = water_available - Qs2;

        S(i) = water_available;
        
        Q(i) = Qs(i) + Qg(i);
    end
end
