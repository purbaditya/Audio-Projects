
%% Smoothing

H_original      = 20*log10(abs(fft(h_original,NFFT)));
H_coarse        = H_original(1:NFFT/2+1);
H_coarse        = H_coarse - mean_coarse;
fvec            = (0:NFFT/2) * fs/NFFT;
if N_smooth > 0
	H_smooth        = smoothSpectrum(H_coarse,fvec,N_smooth);
else
    H_smooth        = H_coarse;
end
H_smooth_label  = [H_smooth fliplr(H_smooth(2:end-1))]/20;

%% Initialization

numpeakfilters  = 10;
numfilters      = numpeakfilters + 2;
global_G        = 1;

% mean shift
s               = h_original*10^(-mean_coarse/20);
n               = numel(s);
K               = 1;
x               = [1 zeros(1,K*n-1)];

% error placeholders
error           = [];
errorepoch      = [];

% generate label
yl              = s;

%% Initialize parameters

paramstruct         = initializeparams_v3(numfilters,K*n,param_approx);
paramstruct.y{1}    = x;    % input
paramstructinit     = initializeparams_v3(numfilters,K*n,param_approx);

yi      = adaptive_cascadefiltforward(fs, paramstruct, numfilters);
X       = 20*log10(abs(fft(yi,NFFT)));
X       = X(1:NFFT/2+1);
YL      = 20*log10(abs(fft(yl,NFFT)));
YL      = YL(1:NFFT/2+1);

% Plot initialization
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

%% calculate error
err     = H_smooth - X;
