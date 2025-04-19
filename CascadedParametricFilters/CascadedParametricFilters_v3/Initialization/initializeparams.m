function paramstruct = initializeparams(numfilters,numsamples,filtertype,fc_initial,G_initial,fb_initial)

l = 2; % order of the peak filter
paramstruct             = {};
paramstruct.filtertype  = cell(numfilters,1);
paramstruct.y           = cell(numfilters+1,1);
paramstruct.ap_y        = cell(numfilters,1);
paramstruct.xh          = cell(numfilters,1);
paramstruct.G           = cell(numfilters,1);
paramstruct.fc          = cell(numfilters,1);
paramstruct.fb          = cell(numfilters,1);
paramstruct.dydx        = cell(numfilters+1,1);
paramstruct.dydG        = cell(numfilters,1);
paramstruct.dydfc       = cell(numfilters,1);
paramstruct.dydfb       = cell(numfilters,1);
paramstruct.dxhdG       = cell(numfilters,1);
paramstruct.dxhdfc      = cell(numfilters,1);
paramstruct.dxhdfb      = cell(numfilters,1);

% for adam optimization
paramstruct.mG      = cell(numfilters,1);
paramstruct.vG  	= cell(numfilters,1);
paramstruct.mfc   	= cell(numfilters,1);
paramstruct.vfc     = cell(numfilters,1);
paramstruct.mfb     = cell(numfilters,1);
paramstruct.vfb     = cell(numfilters,1);

% for diffgrad optimization
paramstruct.gfc         = cell(numfilters,1);
paramstruct.gfb         = cell(numfilters,1);
paramstruct.gG          = cell(numfilters,1);
paramstruct.dydx{end}   = zeros(1,numsamples);

for i = 1 : numfilters
    
    % initial values
    paramstruct.G{i}    = G_initial{i};
    paramstruct.fc{i}   = fc_initial{i};
    paramstruct.fb{i}   = fb_initial{i};
    paramstruct.y{i}    = zeros(1,numsamples);
    
    paramstruct.filtertype{i}   = filtertype{i}; 
    paramstruct.ap_y{i}         = 0;
    paramstruct.dydx{i}         = zeros(1,numsamples);
    paramstruct.dydG{i}         = zeros(1,numsamples);
    paramstruct.dydfc{i}        = zeros(1,numsamples);
    paramstruct.dydfb{i}        = zeros(1,numsamples);
    
    % for adam optimization
    paramstruct.mG{i}   = 0;
    paramstruct.vG{i}   = 0;
    paramstruct.mfc{i}  = 0;
    paramstruct.vfc{i}  = 0;
    paramstruct.mfb{i}  = 0;
    paramstruct.vfb{i}  = 0;
    
    % for diffgrad optimization
    paramstruct.gfc{i}  = 0;
    paramstruct.gfb{i}  = 0;
    paramstruct.gG{i}   = 0;    
    
    paramstruct.xh{i}       = zeros(1,numsamples+l);
    paramstruct.dxhdfc{i}   = zeros(1,numsamples+l);
    paramstruct.dxhdfb{i}   = zeros(1,numsamples+l);
    paramstruct.dxhdx{i}    = zeros(1,numsamples+l);
    paramstruct.dxhdG{i}    = zeros(1,numsamples+l);
    
end

paramstruct.y{end}  = zeros(1,numsamples);