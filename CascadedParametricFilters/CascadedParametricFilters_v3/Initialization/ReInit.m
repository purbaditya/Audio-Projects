% reinitialize peak filter gains
start   = 2;
k1      = 0.5;
k2      = 0.1;

for i = start : numfilters-1
    filtertype{i} = 'peak';
    fc_initial{i} = loc(i);
    if i == 2
        G_initial{i}  = G_initial{i} + k2 * E(i-1);
    else
        if loc(i) > 7000
            G_initial{i}  = G_initial{i} + k2 * E(i-1);
        else
            G_initial{i}  = G_initial{i} + k2 * E(i-1);
        end
    end
    fb_initial{i} = bw(i);
end

disp('done');
numpeakfilters  = numfilters-2;

%% mean shift
s       = h_original*10^(-mean_coarse/20);
n       = numel(s);
K       = 1;
x       = [1 zeros(1,K*n-1)];

% error placeholders
error       = [];
errorepoch  = [];

%% generate label

yl      = s;

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
semilogx(fvec,H_smooth,'Color',hsublack,'LineWidth',1.5);hold on;
semilogx(fvec,X,'Color',hsucyan,'LineWidth',1.5);
grid on;
xlabel('Hz','Interpreter','latex');
ylabel('dB','Interpreter','latex');
legend('Desired','Initial','Interpreter','latex','AutoUpdate','off','location','northwest');
xlim([20 20000]);

for i = 2 : numel(paramstructinit.fc)-1
    plot([paramstructinit.fc{i} paramstructinit.fc{i}],ylim,'Color',rgbval,'LineWidth',0.25);
end
hold off;
title(strcat('No. of Peak Filters : ',num2str(numpeakfilters)));
set(gca,'fontname','times');
set(findall(gcf,'-property','FontSize'),'FontSize',20);

figure(102)
fcn_visualize_approximation(yl',yi',paramstructinit,fs,'');
set(gca,'fontname','times');
set(findall(gcf,'-property','FontSize'),'FontSize',40);
set(gca,'linewidth',2.5);

E_init      = mean((H_smooth(2:end)-X(2:end)).^2);