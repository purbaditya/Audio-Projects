function paramstruct = initializeparamsforwards(numfilters,numsamples,filtertype,fc_desired,G_desired,fb_desired)
l = 2; % order of the peak filter
paramstruct = {};
paramstruct.filtertype = cell(numfilters,1);
paramstruct.y = cell(numfilters+1,1);
paramstruct.ap_y = cell(numfilters,1);
paramstruct.xh = cell(numfilters,1);
paramstruct.G = cell(numfilters,1);
paramstruct.fc = cell(numfilters,1);
paramstruct.fb = cell(numfilters,1);

for i = 1:numfilters
    
    %% desired values %%%%%%%%
    paramstruct.G{i} = G_desired{i};%20-5*i;
    paramstruct.fc{i} = fc_desired{i};%500*i;
    paramstruct.fb{i} = fb_desired{i};%80*i;
    %% %%%%%%%%%%%%%%%%%%%%%%%
    paramstruct.filtertype{i} = filtertype{i}; 
    paramstruct.y{i} = zeros(1,numsamples);
    paramstruct.ap_y{i} = 0;
    paramstruct.xh{i} = zeros(1,numsamples+l);
end
paramstruct.y{end} = zeros(1,numsamples);