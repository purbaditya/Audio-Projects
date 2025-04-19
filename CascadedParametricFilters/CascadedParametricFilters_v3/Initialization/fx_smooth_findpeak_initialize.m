%% Parameters

mpp             = 0.01;                                 % maximum peak prominence threshhold

%% Smoothing

H_original      = 20*log10(abs(fft(h_original,NFFT)));
H_coarse        = H_original(1:NFFT/2+1);
mean_coarse     = mean(H_coarse);
H_coarse        = H_coarse - mean_coarse;
fvec            = (0:NFFT/2) * fs/NFFT;
if N_smooth > 0
	H_smooth        = smoothSpectrum(H_coarse,fvec,N_smooth);
	H_smooth_label  = smoothSpectrum(H_coarse,fvec,N_smooth);
else
    H_smooth        = H_coarse;
    H_smooth_label  = H_coarse;
end
H_smooth_label = [H_smooth_label fliplr(H_smooth_label(2:end-1))]/20;       % Wofür ist das?


%% Find Peaks

% peaks absolute value
[vp,locabs]     = findpeaks(abs(H_smooth),'MinPeakProminence',mpp);

% peaks normal value
[v,locnorm]     = findpeaks(H_smooth,'MinPeakProminence',mpp);

% identify notch without width reference
inv_H_smooth    = max(H_smooth) - H_smooth;
[vn,locminnorm] = findpeaks(inv_H_smooth,'MinPeakProminence',mpp);

% find unique peaks and notches (no repetition)
[~,p]       = find(v<0);        % Peaks in the negative half-plane
p           = locnorm(p);
valn        = H_smooth(locminnorm);
[~,n]       = find(valn>0);     % Notches in the positive half-plane
n           = locminnorm(n);
loc         = sort(unique([locabs locnorm locminnorm]));
pn          = zeros(1,numel(loc));

for i = 1:length(p)
    p(i) = find(loc==p(i));     % Indices of negative peaks in list of peaks/notches
end
for i = 1:length(n)
    n(i) = find(loc==n(i));     % Indices of positive notches in list of peaks/notches
end

% covert location indices to frequencies and store frequencies and magnitudes
val         = H_smooth(loc);
loc1        = fvec(loc);

%% adjust/delete negative peaks and positive notches

% Peaks
[~,idx]         = find(val(p)>=-5);
val(p(idx))     = 1;
[~,idx1]        = find(val(p)<-5);
val(p(idx1))    = [];
loc(p(idx1))    = [];

% Notches
[~,idx]         = find(val(n)<=30);
val(n(idx))     = -1;
[~,idx2]        = find(val(n)>=30);
val(n(idx2))    = [];
loc(n(idx2))    = [];

if loc(end) > 511
    loc(end) = [];
    val(end) = [];
end

%% initialize bandwidths of peak filters

bw      = zeros(size(loc));

for index = 1 : numel(loc)
    k   = 1;
    if H_smooth(loc(index)) > H_smooth(loc(index)-k)    % peak
        L   = 2;
        if abs(H_smooth(loc(index)) - H_smooth(loc(index)-k)) > abs(H_smooth(loc(index)) - H_smooth(loc(index)+k))
            while H_smooth(loc(index)-k) > H_smooth(loc(index)) - L
                k = k+1;
                if (H_smooth(loc(index)-k) > H_smooth(loc(index)-k+1)) || (loc(index) <= k+1)
                    break;
                end
            end
        else
            while H_smooth(loc(index)+k) > H_smooth(loc(index)) - L
                k = k+1;
                if (H_smooth(loc(index)+k) > H_smooth(loc(index)+k-1)) || (loc(index) + k >= length(H_smooth)-2)
                    break;
                end
            end
        end
    else        % notch
        L = 1.5;
        if abs(H_smooth(loc(index)) - H_smooth(loc(index)-k)) > abs(H_smooth(loc(index)) - H_smooth(loc(index)+k))
            while H_smooth(loc(index)-k) < H_smooth(loc(index)) + L
                k = k+1;
                if (H_smooth(loc(index)-k) < H_smooth(loc(index)-k+1)) || (loc(index) <= k+1)
                    break;
                end
            end
        else
            while H_smooth(loc(index)+k) < H_smooth(loc(index)) + L && loc(index)+ k < length(H_smooth)-2
                k = k+1;
                if (H_smooth(loc(index)+k) < H_smooth(loc(index)+k-1)) || (loc(index) + k >= length(H_smooth)-2)
                    break;
                end
            end
        end
    end
    bw(index) = 2*k;
end

bw      = bw * fs/NFFT;
loc2    = loc;
loc     = loc * fs/NFFT;

%% delete/adjust centre freq and bandwidths based on peak proximity

i       = 1;        % counter
f_th    = 300;      % minimum proximity in Hz
M_th    = 2;        % additional magnitude difference threshhold in dB

while i < length(loc) - 1
    if (loc(i+1) < loc(i) + f_th) && val(i)>0 && val(i+1)<0 && (abs(val(i)-H_smooth(round(NFFT*loc(i+1)/fs)))<=M_th)
        loc(i+1)    = [];
        loc2(i+1)   = [];
        val(i+1)    = [];
        bw(i+1)     = [];
        i           = i-1;
    end
    i = i+1;
end

i       = 1;        % counter
while i < length(loc)-1
    if (loc(i+1) < loc(i) + f_th) && val(i)<0 && val(i+1)>0 && (abs(val(i+1)-H_smooth(round(NFFT*loc(i)/fs)))<=M_th)
        loc(i)      = [];
        loc2(i)     = [];
        val(i)      = [];
        bw(i)       = [];
        i           = i-1;
    end
    i = i+1;
end

i       = 1;        % counter
while i < length(loc)-1
    if (loc(i+1) < loc(i) + f_th) && val(i)>0 && val(i+1)>0 && (abs(val(i)-val(i+1))<=M_th)
        loc(i)      = loc(i)+(loc(i+1) - loc(i))*0.5;
        bw(i)       = bw(i) + bw(i+1);
        loc(i+1)    = [];
        loc2(i+1)   = [];
        val(i+1)    = [];
        bw(i+1)     = [];
        i           = i-1;
    end
    i = i+1;
end

%% initialize cut-off frquencies of shelving filters (not very good)

H_smooth2   = H_smooth_label;
binfc       = 2;

% fc LFS
if H_smooth2(1) < H_smooth2(binfc)
    while (20*H_smooth2(binfc) < 20*H_smooth2(1)+2) && (binfc < round(800*NFFT/fs))
        binfc = binfc+1;
    end
else
    while (20*H_smooth2(binfc) > 20*H_smooth2(1)-2) && (binfc < round(800*NFFT/fs))
        binfc = binfc+1;
    end
end
if binfc > round(850*NFFT/fs)
    binfc = 2;
end
locfcls = binfc * fs/NFFT;

binfc = 12;
% fc HFS
if H_smooth2(end-binfc) > H_smooth2(end-1)
    while (20*H_smooth2(end-binfc) < 20*H_smooth2(end-1)+1) && (binfc < round(2150*NFFT/fs))
        binfc = binfc+1;
    end
else
    while (20*H_smooth2(end-binfc) > 20*H_smooth2(end-1)-1) && (binfc < round(2150*NFFT/fs))
        binfc = binfc+1;
    end
end
if binfc > round(2150*2*length(H_smooth)/fs)
    binfc = 12;
end
locfchs = fs/2 - binfc * fs/NFFT;

%% initialize gain

Init;

%% calculate error
err     = H_smooth - X;
E       = err(loc2);

%% re-initialize gain based on error
ReInit;
erra = mean(abs(H_smooth(2:464)-X(2:464)));

while erra > 2 % careful with the value, might run in infinite loop, this can be done as pre or post-processing in a for-loop as well
    % calculate error
    err     = H_smooth(2:end-1)-X(2:end-1);
    E       = err(loc2);
    % re-initialize gain
    ReInit;
    erra    = mean(abs(H_smooth(2:464)-X(2:464)));
end
disp('done');