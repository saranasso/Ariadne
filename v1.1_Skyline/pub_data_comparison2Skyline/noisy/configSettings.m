%% config
global settings
settings.deltaLabelR=10.0083;
settings.deltaLabelK=8.0142;
settings.labelName='AQUA_heavy';
settings.invertedLabel=0; % if the reference is the heavy =0, otherwise =1
settings.constantTarget = 'light'; % it matters only for dilution series (settings.dilution_series=1). If heavy is diluted then constant target is the light, otherwise if light is changing for any reason and heavy is kept constant (like in real quant experiments) then the constant target is the heavy
settings.multFactor = 10; % it corrects for possible ( but discouraged) multiplicative differences in the constantTarget (i.e., reference) amount. If the reference amount in your real sample is N times less than the amount spiked in the calibration experiment then you want to multiply by N.
settings.ionsResolution=0.02; % it depends on the transition list resolution (e.g.,  minimum distance between fragment ions, here was 0.02 Da)
settings.dilution_series=1;
settings.calibration=0;
settings.interferingBackground = 0;
settings.absoluteQuantification=1;
settings.timeUnit='min'; % sec or min are the only 2 allowed values. It refers to the time unit used in the metadata (it must be consistent!)
settings.mProphet=1; % if any mProphet output is used as an input
settings.mProphetFilter=0; % if you want to filter targets by mProphet results
settings.FDRthreshold=0.05; % FDR cutoff for mProphet output
settings.massShiftTolerance=0.001; % it should account for rounding errors: be sure it's smaller than settings.ionsResolution
settings.minXcorr=0.9; % matching probabilty: do not modify it, this is the value to be used for transitions' pair matching
settings.xCorrSignificativity=0.01; %p-value for matching significance
settings.noiseLevel=100; % used as a scaling factor, suggested as 100
settings.minExpectPeakHeight=100;
settings.expPeakWidth=60; % ~twice the chromatography peak width
settings.filterTrans=1; % if settings.mProphetFilter=0 this MUST be 1 otherwise it can be also 0 and all acquired transitions, also not present in he input transition list will be used
settings.interferingBackground = 0; % set it to 1 if and only if you are diluting the heavy peptides and using as a constant reference the light version in a complex background where the same peptide can be present

if settings.constantTarget == 'heavy'
    settings.interferingBackground = 0; 
end

if settings.mProphetFilter == 0 
    settings.filterTrans=1;
end

settings.Skyline=1; % if any mProphet output is used as an input
settings.SkylineFilter=0; % if you want to filter targets by mProphet results