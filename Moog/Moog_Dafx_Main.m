clear all;close all;clc;
% P. Bhattacharya
% Function to test Moog filter parameter adaptation

%% Params
fc      = 1000;                   % initial cut-off frequency (to be adapted or learned)
fc_l    = 3500;                     % ground truth or label cut-off frequency
k       = 0.2;                     % feedback coeff (Q-factor?) - to be adapted as well, not yet implemented
k_l     = 0.5;
fs      = 44100;                   % sampling rate
T       = 1;                       % time period in sec
N       = 4;                       % size of Moog ladder
%x       = zeros(1,T*fs); x(1) = 1; % input dirac
x       = 0.5*(2*rand(1,T*fs)-1);  % input dirac
y       = zeros(1,length(x));      % output
e       = ones(1,length(x));      % error
osf     = 2;                       % oversampling factor
eta_p   = 0.02;                    % learning rate
eta_k   = 0.001;                    % learning rate

% Oversampling
x2    = zeros(1,osf*length(x)); x2(1:osf:end) = x;
p     = 2*pi*fc/(fs*osf);
p_l   = 2*pi*fc_l/(fs*osf);
y2    = zeros(1,osf*length(y));
y2_l  = zeros(1,osf*length(y));
dy2dp = zeros(1,osf*length(y));
dy2dk = zeros(1,osf*length(y));

% Anti-imaging filter
h = fir1(10,0.5);
x2 = filter(h,1,x2);

%% Initialization and ground-truth generation
% States
xd_old  = 0;
y_old   = zeros(N,1);
y_temp  = zeros(N,1);
dxd_old_dp = 0;
dy_old_dp  = zeros(N,1);
dy_temp_dp = zeros(N,1);
for n = 1:length(x2)
    [y_temp,y_old,xd_old] = Moog_Dafx_Arch(x2(n),y_temp,y_old,xd_old,k_l,p_l,N);
    % Label or ground truth signal
    y2_l(n) = y_temp(end);               
end

%% Training
p_vec = [];k_vec = [];epoch = 1;

% Stop criterion
while epoch < 2000 && abs(mean(e(1:fs/10))) > 4e-8% Tolerance - 0.001
    xd_old  = 0;
    y_old   = zeros(N,1);
    y_temp  = zeros(N,1);
    dxd_old_dp = 0;
    dy_old_dp  = zeros(N,1);
    dy_temp_dp = zeros(N,1);
    dxd_old_dk = 0;
    dy_old_dk  = zeros(N,1);
    dy_temp_dk = zeros(N,1);
    % *No need to iterate across all samples, speeds up convergence
    for n = 1:length(x2)/10             
        [y_temp,dy_temp_dp,dy_temp_dk,y_old,dy_old_dp,dy_old_dk,xd_old,dxd_old_dp,dxd_old_dk] = Moog_Dafx_Arch_forward_backward(x2(n),y_temp,dy_temp_dp,dy_temp_dk,y_old,dy_old_dp,dy_old_dk,xd_old,dxd_old_dp,dxd_old_dk,k,p,N);
        y2(n) = y_temp(end);
        
        % Error gradient
        e(n) = (abs(y2(n))-abs(y2_l(n)))*sign(y2(n));
        
        % Collect all gradients,
        % Trick - factor (fs/n^4) skews the error, only initial samples are responsible for gradient
        % avoids occasional wrong gradients
        if epoch < 250
            dy2dp(n) = dy_temp_dp(end)*e(n)*(1.1*fs/n^5);
            dy2dk(n) = dy_temp_dk(end)*e(n);
        elseif epoch < 400
            dy2dp(n) = dy_temp_dp(end)*e(n)*(0.5*fs/n^4);
            dy2dk(n) = dy_temp_dk(end)*e(n);
        elseif epoch < 1200
            dy2dp(n) = dy_temp_dp(end)*e(n)*(fs/n^3);
            dy2dk(n) = dy_temp_dk(end)*e(n);
            eta_p = 0.001;
            eta_k = 0.001;
        else
            dy2dp(n) = dy_temp_dp(end)*e(n)*(1.1*fs/n^3);
            dy2dk(n) = dy_temp_dk(end)*e(n);
            eta_p = 0.0005;
            eta_k = 0.0005;
        end
        
        % Instanteneous update (backpropagation)
        %p = p-eta*dy2dp(n);
    end
    % Visualize update per epoch
    figure(9);plot(y2,'r');hold on;plot(y2_l,'k');grid;xlim([0 500]);ylim([-0.2 0.2]);
    xlabel('n');ylabel('y(n)');legend('prediction','reference');hold off;
    
    % Normal update with gradient clipping
    gradval_p = min(0.25,abs(sum(dy2dp(:))))*sign(sum(dy2dp(:)));   
    gradval_p = max(0.0005,abs(gradval_p))*sign(gradval_p);
    if epoch < 500
        gradval_k = min(1,abs(sum(dy2dk(:))))*sign(sum(dy2dp(:)));
    else
        gradval_k = sum(dy2dk(:));
    end
    % update
    p = p-eta_p*gradval_p;
    k = k-eta_k*gradval_k;
    p_vec = [p_vec p];
    k_vec = [k_vec k];
    epoch = epoch+1;
end

% fc = p*fs*osf/(2*pi);
%% Plots
figure(10);
grid on;
plot(p_vec*fs*osf/(2*pi),'r');hold on;plot(p_l*fs*osf/(2*pi)*ones(1,length(p_vec)),'k');hold off;grid on;
xlabel('Epochs');ylabel('fc in Hz');xlim([1 epoch]);legend('prediction','reference');

figure(11);
grid on;
plot(k_vec,'r');hold on;plot(k_l*ones(1,length(p_vec)),'k');hold off;grid on;
xlabel('Epochs');ylabel('k');xlim([1 epoch]);legend('prediction','reference');

% Antialiasing filter at fs2
y2 = filter(h,1,y2);

y = osf*y2(1:osf:end);
% figure(1);plot(y);

fc_new = p*fs*osf/(2*pi);
[h,w] = freqz(y,1,8192,fs);
figure(2);plot(w,20*log10(abs(h)));hold on;
plot([fc_new, fc_new], [-40, 30]);
hold off;
ax = gca;
ax.YLim = [-40 20];
ax.XLim = [20 fs/2];
ax.XScale = 'log';
xlabel('Frequency in Hz');
ylabel('Magnitude (dB)');
grid on;