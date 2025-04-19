function yl = adaptive_cascadefiltforward(fs, paramstruct, numlayers)

% read from the structure
filtertype      = paramstruct.filtertype;
y               = paramstruct.y;
xh              = paramstruct.xh;
G               = paramstruct.G;
fc              = paramstruct.fc;
fb              = paramstruct.fb;

for n = 1:length(y{1})
    % forward calculation per sample
    for m = 1 : 1 : numlayers
        [y{m+1}, xh{m}, ~] = filtlayer(y{m},y{m+1},xh{m},fc{m},fb{m},G{m},filtertype{m},n,fs);
    end
end

% copy back to structure
paramstruct.filtertype  = filtertype;
paramstruct.y           = y;
paramstruct.xh          = xh;
paramstruct.G           = G;
paramstruct.fc          = fc;
paramstruct.fb          = fb;

yl 	= y{end};

end