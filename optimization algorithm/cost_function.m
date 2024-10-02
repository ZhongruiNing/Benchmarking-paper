function CF = cost_function(x, param, observation, fn_hm)
    addpath(genpath("../"))
    simulation = fn_hm(x, param);

    NSE = nash_sutcliffe_efficiency(observation, simulation);
    RE = abs(relative_error(observation, simulation));
    % KGE = klinggupta(observation, simulation);

    CF = 1 - NSE + abs(RE);
end
