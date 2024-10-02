function KGE = klinggupta(observation, simulation)
    %Nash Sutcliffe Efficiency measure
    simulation(isnan(observation)) = NaN;
    observation(isnan(simulation)) = NaN;

    cflow = [simulation, observation];

    sdsimulation = nanstd(simulation);
    sdobservation = nanstd(observation);
    
    msimulation = nanmean(simulation);
    mobservation = nanmean(observation);

    r = corrcoef(cflow, 'rows', 'pairwise'); 
    r = r(1, 2);
    relvar = sdsimulation / sdobservation;
    bias = msimulation / mobservation;

    %KGE timeseries 
    KGE = 1 - sqrt(((r - 1) ^ 2) + ((relvar - 1) ^ 2) + ((bias - 1) ^ 2)); 
end