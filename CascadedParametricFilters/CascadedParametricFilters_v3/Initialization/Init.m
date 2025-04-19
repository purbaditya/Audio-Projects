%% Initialize

k           = 0.4;           % initial gain scale
numfilters  = numel(loc);    % number of peak filters

% LFS
filtertype  = {'lowshelving'};
fc_initial  = {locfcls};    
G_initial   = {H_smooth(1)};
fb_initial  = {0};
loc         = [locfcls loc];
val         = [H_smooth(1) val];
bw          = [0 bw];
numfilters  = numfilters + 1;

% Peak filters
start       = 2;
for i = start : numfilters
    filtertype{i}   = 'peak';
    fc_initial{i}   = loc(i);
    if i == 2           % All gains are similarily calculated at the moment
        G_initial{i}  = k * val(i);
    else
        if loc(i) > 10000
            G_initial{i}  = k * val(i);
        else
            G_initial{i}  = k * val(i);
        end
    end
    fb_initial{i}   = bw(i);
end

% HFS
numfilters      = numfilters+1;
i               = numfilters;
filtertype{i}   = 'highshelving';
fc_initial{i}   = locfchs;
G_initial{i}    = H_smooth(end-40);
fb_initial{i}   = 0;

global_G        = 1; % not used currently

numpeakfilters  = numfilters - 2;

%% mean shift

s               = h_original*10^(-mean_coarse/20);
n               = numel(s);
K               = 1;
x               = [1 zeros(1,K*n-1)];

% error placeholders
error           = [];
errorepoch      = [];

%% generate label

yl              = s;

%% initialize network parameters

paramstruct         = initializeparams(numfilters,K*n,filtertype,fc_initial,G_initial,fb_initial);
paramstruct.y{1}    = x;    % input
paramstructinit     = initializeparams(numfilters,K*n,filtertype,fc_initial,G_initial,fb_initial);

%% plot initial HRTF

yi      = adaptive_cascadefiltforward(fs, paramstruct, numfilters);
X       = 20*log10(abs(fft(yi,NFFT)));
X       = X(1:NFFT/2+1);
YL      = 20*log10(abs(fft(yl,NFFT)));
YL      = YL(1:NFFT/2+1);

figure(2);
semilogx(fvec,X); hold on;
semilogx(fvec,H_smooth); 
semilogx(fvec,H_smooth-X);
grid on;
xlabel('Frequency in Hz','Interpreter','latex');
ylabel('Magnitude in dB','Interpreter','latex');
legend('initial','desired','error','Interpreter','latex','AutoUpdate','off');
axis([20 20000 -30 20]);
for i = 2:numel(paramstructinit.fc)-1
    plot([paramstructinit.fc{i} paramstructinit.fc{i}],ylim,'Color',rgbval,'LineWidth',1.5);
end
hold off;
title(strcat('No. of Peak Filters : ',num2str(numpeakfilters)));
