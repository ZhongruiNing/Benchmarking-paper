function [Q, E, SN] = YWBMbS(x, params)

    P   = x(:, 1);
    T   = x(:, 2);
    PET = x(:, 3);
    
    Ks      = params(1);
    Kg      = params(2);
    beta    = params(3);
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
    
    Pr      = zeros(n, 1);
    Ps      = zeros(n, 1);
    
    % Main body for WBM
    for i = 2 : n
        S0  = S(i - 1);
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
        check_complex(Mt)

        SN(i) = SN0 + SFt - Mt;
        P(i) = RFt + Mt;

        if T(i) < Tsnow
            water_available = S0;

            Qg(i) = 0;
            Qs(i) = P(i);

            Q(i) = Qg(i) + Qs(i);

            temp_x = PET(i) / 0.1;
            w = 1 / (1 - alpha);
            actual_evap = water_available * (1 + temp_x - (1 + temp_x ^ w) ^ (1 / w));
            E(i) = E(i) + actual_evap;

            water_available = max(0.001, water_available - actual_evap);
            
            S(i) = water_available;
        else
            % Underground runoff based on runoff in last time step
            Qg(i) = Kg * S0;
            water_available = S0 - Qg(i) + P(i);
            
            % Evaporation Module
            temp_x = PET(i) / water_available;
            w = 1 / (1 - alpha);
            actual_evap = water_available * (1 + temp_x - (1 + temp_x ^ w) ^ (1 / w));
            E(i) = E(i) + actual_evap;

            water_available = max(0.001, water_available - actual_evap);

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
end
