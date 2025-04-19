function [dydx, dydG, dydfc, dydfb, dxhdx, dxhdG, dxhdfc, dxhdfb] = filtlayerbackprop(x,ap_y,xh,fc,fb,G,dydx, dydG, dydfc, dydfb,dxhdx,dxhdG,dxhdfc,dxhdfb,delta,filtertype,n,fs)
%% Gain, fc, fb, x Adaptation %%
% ------------------------------

K_G_p = 1; % gradient multiplier
K_G_n = 2; % gradient multiplier
K_fb = 1; % gradient multiplier
K_fc = 1; % gradient multiplier
K_x = 1; % gradient multiplier

switch(filtertype)
    case 'peak'
        l = 2; % filter order
        
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
        if G >= 0
            % grad w.r.t G %
            dxhdG(n+l) = 0;
            dydG(n) = K_G_p*delta(n)*(x(n)-ap_y)*log(10)*V0*(1/40);
            
            % grad w.r.t fb , l = 2 represents the order of the peak filter %
            dcdfb = (2*pi/fs)*((sec(pi*fb/fs)^2)/((tan(pi*fb/fs)+1)^2));%(2*pi/fs)*(sec(2*pi*fb/fs)*(sec(2*pi*fb/fs)-tan(2*pi*fb/fs))); % c = tan(2*pi*fb/fs)-sec(2*pi*fb/fs)
            dxhdfb(n+l) = (xh(n+l-2)+d*xh(n+l-1))*dcdfb - d*(1-c)*dxhdfb(n+l-1) + c*dxhdfb(n+l-2);
            tmp = (-c*dxhdfb(n+l) - (xh(n+l)+d*xh(n+l-1))*dcdfb + d*(1-c)*dxhdfb(n+l-1)+dxhdfb(n+l-2));
            dydfb(n) = K_fb*delta(n)*(-H0/2)*tmp;
        else
            % grad w.r.t G %
            dcdG = -((1/10)*log(10)*V0*tan(pi*fb/fs))/((tan(pi*fb/fs)+V0)^2);% corrected
            dxhdG(n+l) = -d*(1-c)*dxhdG(n+l-1) + (d*xh(n+l-1)+xh(n+l-2))*dcdG + c*dxhdG(n+l-2);
            tmp = -dcdG*(xh(n+l)+d*xh(n+l-1)) - c*dxhdG(n+l) + d*(1-c)*dxhdG(n+l-1) + dxhdG(n+l-2);
            dydG(n) = K_G_n*delta(n)*(((x(n)-ap_y)*log(10)*V0*(1/40))-(0.5*H0*tmp));
            
            % grad w.r.t fb , l = 2 represents the order of the peak filter %
            dcdfb = ((2*pi/fs)*V0*sec(pi*fb/fs)^2)/(tan(pi*fb/fs)+V0)^2; % c = tan(2*pi*fb/fs)-sec(2*pi*fb/fs)
            dxhdfb(n+l) = (xh(n+l-2)+d*xh(n+l-1))*dcdfb - d*(1-c)*dxhdfb(n+l-1) + c*dxhdfb(n+l-2);
            tmp = (-c*dxhdfb(n+l) - (xh(n+l)+d*xh(n+l-1))*dcdfb + d*(1-c)*dxhdfb(n+l-1)+ dxhdfb(n+l-2));
            dydfb(n) = K_fb*delta(n)*(-H0/2)*tmp;
        end
        
        % grad w.r.t fc , l = 2 represents the order of the peak filter %
        dddfc = sin(2*pi*fc/fs)*(2*pi/fs);
        dxhdfc(n+l) = -(1-c)*xh(n+l-1)*dddfc - d*(1-c)*dxhdfc(n+l-1) + c*dxhdfc(n+l-2);
        tmp = (-c*dxhdfc(n+l) + (1-c)*xh(n+l-1)*dddfc + d*(1-c)*dxhdfc(n+l-1)+ dxhdfc(n+l-2));
        dydfc(n) = K_fc*delta(n)*(-H0/2)* tmp;
        
        % grad w.r.t x  , l = 2 represents the order of the peak filter %
        %dxhdx(n+l) = 1 - d*(1-c)*dxhdx(n+l-1) + c*dxhdx(n+l-2);
        %dydx(n) = K_x*delta(n)*(1+H0/2 - H0/2 *(-c*dxhdx(n+l) + d*(1-c)*dxhdx(n+l-1)+dxhdx(n+l-2)));
        dydx(n) = K_x*delta(n)*(1+H0/2*(1 + c)); %corrected
    case 'lowshelving'
        l = 2;
        
        Wc = 2 * pi * fc / fs;
        V0 = 10^(G/20);
        H0 = V0 - 1;
        if G >= 0
            c = (tan(Wc/2)-1) / (tan(Wc/2)+1);     % boost
        else
            c = (tan(Wc/2)-V0) / (tan(Wc/2)+V0);   % cut
        end
        
        d = 0;
        
        if G >= 0
            % grad w.r.t G %
            dxhdG(n+l) = 0;
            dydG(n) = K_G_p*delta(n)*0.5*(x(n)+ap_y)*log(10)*V0*(1/20);
            
            % grad w.r.t fc , l = 1 represents the order of the shelving filter %
            dcdfc = (2*pi/fs)*(sec(pi*fc/fs)^2)/((tan(pi*fc/fs)+1)^2);%(2*pi/fs)*(sec(2*pi*fc/fs)*(sec(2*pi*fc/fs)-tan(2*pi*fc/fs))); % c = tan(2*pi*fb/fs)-sec(2*pi*fb/fs)
            dxhdfc(n+l) = -c*dxhdfc(n+l-1) - xh(n+l-1)*dcdfc;
            dydfc(n) = K_fb*delta(n)*(H0/2)*(c*dxhdfc(n+l) + xh(n+l)*dcdfc + dxhdfc(n+l-1));
            
            dydfb(n) = 0;
        else
            % grad w.r.t G %
            dcdG = -((1/10)*log(10)*V0*tan(pi*fc/fs))/((tan(pi*fc/fs)+V0)^2);
            dxhdG(n+l) = -c*dxhdG(n+l-1) - xh(n+l-1)*dcdG;
            tmp = c*dxhdG(n+l) + xh(n+l)*dcdG + dxhdG(n+l-1);
            dydG(n) = K_G_n*delta(n)*(((x(n)+ap_y)*log(10)*V0*(1/40)) + ((H0/2)*tmp));
            
            % grad w.r.t fc , l = 1 represents the order of the shelving filter %
            dcdfc = ((2*pi/fs)*V0*(sec(pi*fc/fs)^2))/((tan(pi*fc/fs)+V0)^2);
            dxhdfc(n+l) = -c*dxhdfc(n+l-1) - xh(n+l-1)*dcdfc;
            dydfc(n) = K_fb*delta(n)*(H0/2)*(c*dxhdfc(n+l) + xh(n+l)*dcdfc + dxhdfc(n+l-1));
            
            dydfb(n) = 0;
        end
        
        
        % grad w.r.t x  , l = 1 represents the order of the shelving filter %
        %dxhdx(n+l) = 1 - c*dxhdx(n+l-1);
        %dydx(n) = K_fc*delta(n)*(1+H0/2 + H0/2 *(c*dxhdx(n+l) + dxhdx(n+l-1)));
        dydx(n) = K_fc*delta(n)*(1+H0/2 + c*H0/2); %corrected
    case 'highshelving'
        l = 2;
        
        Wc = 2 * pi * fc / fs;
        V0 = 10^(G/20);
        H0 = V0 - 1;
        if G >= 0
            c = (tan(Wc/2)-1) / (tan(Wc/2)+1);     % boost
        else
            c = (V0*tan(Wc/2)-1) / (V0*tan(Wc/2)+1);     %(tan(Wb/2)-V0) / (tan(Wb/2)+V0);   % cut
        end
        
        d = 0;
        
        if G >= 0
            % grad w.r.t G %
            dxhdG(n+l) = 0;
            dydG(n) = K_G_p*delta(n)*0.5*(x(n)-ap_y)*log(10)*V0*(1/20);
            
            % grad w.r.t fc , l = 1 represents the order of the shelving filter %
            dcdfc = (2*pi/fs)*((sec(pi*fc/fs)^2)/((tan(pi*fc/fs)+1)^2));%(2*pi/fs)*(sec(2*pi*fc/fs)*(sec(2*pi*fc/fs)-tan(2*pi*fc/fs))); % c = tan(2*pi*fb/fs)-sec(2*pi*fb/fs)
            dxhdfc(n+l) = -c*dxhdfc(n+l-1) - xh(n+l-1)*dcdfc;
            dydfc(n) = K_fb*delta(n)*(-H0/2)*(c*dxhdfc(n+l) + xh(n+l)*dcdfc + dxhdfc(n+l-1));
            
            dydfb(n) = 0;
        else
            % grad w.r.t G %
            dcdG = ((1/10)*log(10)*V0*tan(pi*fc/fs))/((V0*tan(pi*fc/fs)+1)^2);
            dxhdG(n+l) = -c*dxhdG(n+l-1) - xh(n+l-1)*dcdG;
            tmp = c*dxhdG(n+l) + xh(n+l)*dcdG + dxhdG(n+l-1);
            dydG(n) = K_G_p*delta(n)*(((x(n)-ap_y)*log(10)*V0*(1/40)) - ((H0/2)*tmp));
            
            % grad w.r.t fc , l = 1 represents the order of the shelving filter %
            dcdfc = ((2*pi/fs)*V0*sec(pi*fc/fs)^2)/(V0*tan(pi*fc/fs)+1)^2;
            dxhdfc(n+l) = -c*dxhdfc(n+l-1) - xh(n+l-1)*dcdfc;
            dydfc(n) = K_fb*delta(n)*(-H0/2)*(c*dxhdfc(n+l) + xh(n+l)*dcdfc + dxhdfc(n+l-1));
        
            dydfb(n) = 0;
        end
        
        % grad w.r.t x  , l = 1 represents the order of the shelving filter %
        %dxhdx(n+l) = 1 - c*dxhdx(n+l-1);
        %dydx(n) = K_x*delta(n)*(1+H0/2 - H0/2 *(c*dxhdx(n+l) + dxhdx(n+l-1)));
        dydx(n) = K_x*delta(n)*(1+H0/2 - c*H0/2); %corrected
    otherwise
        disp('no such filter exists');
end
% disp(sum(dydx.^2));
% disp(sum(dydG.^2));
% disp(sum(dydfc.^2));
% disp(sum(dydfb.^2));
end
