function [E,e,E_plot] = calculate_error_derivative(y,yi,yl,n,domain,NFFT,H_smooth_label,currentepoch,totalepochs)

    y_temp = [y{end}(1:n) yi(n+1:end)];

    if strcmp(domain,'time')
        % Time Domain
        e       = y_temp-yl;
        E       = e;
        E_plot  = E;
    elseif strcmp(domain,'magnitude')
        % Frequency Domain
        Y       = fft(y_temp,NFFT);
        YL      = fft(yl,NFFT);
        Y_Real  = (real(Y));
        YL_Real = (real(YL));
        Y_Imag  = (imag(Y));
        YL_Imag = (real(YL));

        Y_abs   = (abs(Y).^2);
        YL_abs  = 2*H_smooth_label;         % Warum mal 2?

        if currentepoch < round(0.6*totalepochs)
            E = sign(log10(Y_abs)-YL_abs);  % L1 loss
        else
            E = sign(log10(Y_abs)-YL_abs);  % or L2 loss
        end
        K                           = E(NFFT/2-48+1:NFFT/2+48+1);   
        E(NFFT/2-48+1:NFFT/2+48+1) 	= 0*K;          % or 0.01*K;
        E_plot                      = (log10(Y_abs)-YL_abs);
        E_Real                      = Y_Real.*E./(Y_abs+1e-32);
        E_Imag                      = Y_Imag.*E./(Y_abs+1e-32);
        e                           = real(ifft(E_Real+1i*E_Imag))';

        figure(11);
        semilogx(YL_abs(1:round(length(YL_abs)/2))); grid;
        figure(11);
        plot(Y_Real.^2+Y_Imag.^2); hold on;
        plot(YL_Real.^2+YL_Imag.^2); hold off;
        legend('predicted','desired');
    end
    
end