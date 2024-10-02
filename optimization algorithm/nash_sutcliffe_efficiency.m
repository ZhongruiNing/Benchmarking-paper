function NSE = nash_sutcliffe_efficiency(observation, simulation)
    loc = (~isnan(observation)) & (~isnan(simulation));
    observation = observation(loc);
    simulation = simulation(loc);
    NSE = 1 - sum((observation - simulation) .^ 2) / sum((observation - mean(observation)) .^ 2);
end