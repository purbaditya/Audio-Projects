function [y,xh,ap_y] = filtlayer(x,y,xh,fc,fb,G,filtertype,n,fs)
switch(filtertype)
    case 'peak'
        l = 2;
        Wc = 2 * pi * fc / fs;
        Wb = 2 * pi * fb / fs; % Wb = 2*pi*fb/fs --> fb = Wb*fs/2 = fc/Q
        V0 = 10^(G/20);
        H0 = V0 - 1;
        if G >= 0
            c = (tan(Wb/2)-1) / (tan(Wb/2)+1);     % boost
        else
            c = (tan(Wb/2)-V0) / (tan(Wb/2)+V0);   % cut
        end
        
        d = -cos(Wc);
        
        xh(n+l) = x(n) - d*(1-c)*xh(n+l-1) + c*xh(n+l-2);
        ap_y = -c * xh(n+l) + d*(1-c)*xh(n+l-1) + xh(n+l-2);
        y(n) = 0.5 * H0 * (x(n) - ap_y) + x(n);
    case 'lowshelving'
        l = 1;
        Wc = 2 * pi * fc / fs;
        V0 = 10^(G/20);
        H0 = V0 - 1;
        if G >= 0
            c = (tan(Wc/2)-1) / (tan(Wc/2)+1);     % boost
        else
            c = (tan(Wc/2)-V0) / (tan(Wc/2)+V0);   % cut
        end
        d = 0;
        
        xh(n+l) = x(n) - c*xh(n+l-1);
        ap_y = c * xh(n+l) + xh(n+l-1);
        y(n) = 0.5 * H0 * (x(n) + ap_y) + x(n);
    case 'highshelving'
        l = 1;
        Wc = 2 * pi * fc / fs;
        V0 = 10^(G/20);
        H0 = V0 - 1;
        if G >= 0
            c = (tan(Wc/2)-1) / (tan(Wc/2)+1);     % boost
        else
            c = (V0*tan(Wc/2)-1) / (V0*tan(Wc/2)+1);   % cut
        end
        d = 0;
        
        xh(n+l) = x(n) - c*xh(n+l-1);
        ap_y = c * xh(n+l) + xh(n+l-1);
        y(n) = 0.5 * H0 * (x(n) - ap_y) + x(n);
    otherwise
        disp('no such filter exists');
end
end