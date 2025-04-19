function [s1,s4,m1,m2] = fcn_visualize(paraminit,paramfinal,NFFT,fs,fig_title)

hsured = [0.77255 0 0.25882];
hsugray = [0.5647 0.5176 0.4627];
hsucyan = [0,0.4431,0.4157];

fc_vecfinal      = paramfinal.fc;
G_vecfinal       = paramfinal.G;
fb_vecfinal      = paramfinal.fb;

fc_vecinit      = paraminit.fc;
G_vecinit       = paraminit.G;
fb_vecinit      = paraminit.fb;

p_vec       = numel(fc_vecfinal)-2;
fvec        = ((0:NFFT-1) * fs/NFFT)';

c_gray      = 0.5 * [1, 1, 1];      % Colour for individual peak filters

d           = zeros(NFFT,1);
d(1)        = 1;

m1 = max([G_vecfinal{:} G_vecinit{:}]);
m2 = min([G_vecfinal{:} G_vecinit{:}]);
%% FINAL %%
% Low-Shelving
d_plot      = zeros(NFFT,1);
xh          = zeros(NFFT,1);
for n = 1 : NFFT
    [d_plot,xh,~]   = filtlayer(d,d_plot,xh,fc_vecfinal{1},fb_vecfinal{1},G_vecfinal{1},'lowshelving',n,fs);
end
D_Plot      = fft(d_plot);
D_Plot      = 20*log10(abs(D_Plot));
s1 = semilogx(fvec,D_Plot,'Color',hsured,'LineWidth',2.5); hold on;

d_plot      = zeros(NFFT,1);
xh          = zeros(NFFT,1);
for n = 1 : NFFT
    [d_plot,xh,~]   = filtlayer(d,d_plot,xh,fc_vecfinal{end},fb_vecfinal{end},G_vecfinal{end},'highshelving',n,fs);
end
D_Plot      = fft(d_plot);
D_Plot      = 20*log10(abs(D_Plot));
s2 = semilogx(fvec,D_Plot,'Color',hsured,'LineWidth',2.5); hold on;

for p = 1 : p_vec
    d_plot      = zeros(NFFT,1);
    xh          = zeros(NFFT,1);
    for n = 1 : NFFT
        [d_plot,xh,~]   = filtlayer(d,d_plot,xh,fc_vecfinal{p+1},fb_vecfinal{p+1},G_vecfinal{p+1},'peak',n,fs);
    end
    D_Plot      = fft(d_plot);
    D_Plot      = 20*log10(abs(D_Plot));
    s3 = semilogx(fvec,D_Plot,'Color',hsured,'LineWidth',2.5);
end

%% INIT %%
% Low-Shelving
d_plot      = zeros(NFFT,1);
xh          = zeros(NFFT,1);
for n = 1 : NFFT
    [d_plot,xh,~]   = filtlayer(d,d_plot,xh,fc_vecinit{1},fb_vecinit{1},G_vecinit{1},'lowshelving',n,fs);
end
D_Plot      = fft(d_plot);
D_Plot      = 20*log10(abs(D_Plot));
s4 = semilogx(fvec,D_Plot,'Color',hsucyan,'LineWidth',2.5); hold on;

d_plot      = zeros(NFFT,1);
xh          = zeros(NFFT,1);
for n = 1 : NFFT
    [d_plot,xh,~]   = filtlayer(d,d_plot,xh,fc_vecinit{end},fb_vecinit{end},G_vecinit{end},'highshelving',n,fs);
end
D_Plot      = fft(d_plot);
D_Plot      = 20*log10(abs(D_Plot));
s5 = semilogx(fvec,D_Plot,'Color',hsucyan,'LineWidth',2.5); hold on;

for p = 1 : p_vec
    d_plot      = zeros(NFFT,1);
    xh          = zeros(NFFT,1);
    for n = 1 : NFFT
        [d_plot,xh,~]   = filtlayer(d,d_plot,xh,fc_vecinit{p+1},fb_vecinit{p+1},G_vecinit{p+1},'peak',n,fs);
    end
    D_Plot      = fft(d_plot);
    D_Plot      = 20*log10(abs(D_Plot));
    s6 = semilogx(fvec,D_Plot,'Color',hsucyan,'LineWidth',2.5);
end
end