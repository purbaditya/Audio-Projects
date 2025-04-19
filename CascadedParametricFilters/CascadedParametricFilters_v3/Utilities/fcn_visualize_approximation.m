function [] = fcn_visualize_approximation(h_original,h_approx,param,fs,fig_title)

    hsublack    = [0 0 0];
    hsucyan     = [0,0.4431,0.4157];

    fc_vec      = param.fc;
    G_vec       = param.G;
    fb_vec      = param.fb;

    p_vec       = numel(fc_vec)-2;

    NFFT        = length(h_original);
    fvec        = ((0:NFFT-1) * fs/NFFT)';

    H_Original 	= fft(h_original,NFFT);
    H_Original  = 20*log10(abs(H_Original));
    H_Smooth    = smoothSpectrum(H_Original,fvec,12);

    H_Approx 	= fft(h_approx,NFFT);
    H_Approx    = 20*log10(abs(H_Approx));

    c_gray      = 0.5 * [1, 1, 1];      % Colour for individual peak filters

    d           = zeros(NFFT,1);
    d(1)        = 1;

    % Low-Frequency Shelving
    d_plot      = zeros(NFFT,1);
    xh          = zeros(NFFT,1);
    for n = 1 : NFFT
        [d_plot,xh,~]   = filtlayer(d,d_plot,xh,fc_vec{1},fb_vec{1},G_vec{1},'lowshelving',n,fs);
    end
    D_Plot      = fft(d_plot,NFFT);
    D_Plot      = 20*log10(abs(D_Plot));
    s5          = semilogx(fvec,D_Plot,'Color',c_gray,'LineStyle',':','LineWidth',2); hold on;

    % High-Frequency Shelving
    d_plot      = zeros(NFFT,1);
    xh          = zeros(NFFT,1);
    for n = 1 : NFFT
        [d_plot,xh,~]   = filtlayer(d,d_plot,xh,fc_vec{end},fb_vec{end},G_vec{end},'highshelving',n,fs);
    end
    D_Plot      = fft(d_plot,NFFT);
    D_Plot      = 20*log10(abs(D_Plot));
    semilogx(fvec,D_Plot,'Color',c_gray,'LineStyle',':','LineWidth',2);

    % Peak Filter
    for p = 1 : p_vec
        d_plot      = zeros(NFFT,1);
        xh          = zeros(NFFT,1);
        for n = 1 : NFFT
            [d_plot,xh,~]   = filtlayer(d,d_plot,xh,fc_vec{p+1},fb_vec{p+1},G_vec{p+1},'peak',n,fs);
        end
        D_Plot      = fft(d_plot,NFFT);
        D_Plot      = 20*log10(abs(D_Plot));
        s4 = semilogx(fvec,D_Plot,'Color',c_gray,'LineWidth',2);
    end

    % Approximation
    s2 = semilogx(fvec,H_Smooth,'Color',hsublack,'LineWidth',3);
    s3 = semilogx(fvec,H_Approx,'Color',hsucyan,'LineWidth',3);
    hold off;
    xlim([50 20e3]);
    xlabel('Frequency in Hz','Interpreter','latex'); ylabel('Magnitude in dB','Interpreter','latex'); grid on;
    title(fig_title);
    legend([s2,s3,s4,s5],'Target','Initial','Peak','Shelving','Interpreter','latex','location','northwest');

end