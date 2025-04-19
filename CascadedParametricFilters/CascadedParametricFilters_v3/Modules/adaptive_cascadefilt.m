function [y,paramstruct,global_G,e_plot] = adaptive_cascadefilt(yl, yi, fs, paramstruct, numlayers, lr_multiplier,global_G,e_plot,f_update,updaterule,NFFT,H_smooth_label,optimizeall,currentepoch,totalepochs)

% copy utility
struct2read;

for n = 1 : length(y{1})
    % forward calculation per sample
    for m = 1:1:numlayers
        [y{m+1},xh{m},ap_y{m}] = filtlayer(y{m},y{m+1},xh{m},fc{m},fb{m},G{m},filtertype{m},n,fs);
    end
    
    % error calculation per sample
    [~,e,E_plot]    = calculate_error_derivative(y,yi,yl,n,'magnitude',NFFT,H_smooth_label,currentepoch,totalepochs);
    e_temp          = sum(E_plot.^2)/numel(E_plot);%(y{end}-yl);
    e_plot          = [e_plot e_temp];
    
    % backward calculation per sample
    dydx{end} = e;
    for m = numlayers:-1:1
        [dydx{m}, dydG{m}, dydfc{m}, dydfb{m}, dxhdx{m}, dxhdG{m}, dxhdfc{m}, dxhdfb{m}]...
            = filtlayerbackprop(y{m},ap_y{m},xh{m},fc{m},fb{m},G{m},dydx{m},dydG{m},dydfc{m},dydfb{m},dxhdx{m},dxhdG{m},dxhdfc{m},dxhdfb{m},dydx{m+1},filtertype{m},n,fs);
    end
    
    %% parameter update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for m = 1:1:numlayers
        
        %% G -> peak_filter_gradient_scale=1, shelving_filter_gradient_scale = 0.05 %%
        if m~=numlayers || m~= 1 
            [G{m},mG{m},vG{m},gG{m}]  = optimize(G{m},dydG{m},mG{m},vG{m},gG{m},n,lr_multiplier,1,m,f_update,updaterule,optimizeall);
        else                     
            [G{m},mG{m},vG{m},gG{m}]  = optimize(G{m},dydG{m},mG{m},vG{m},gG{m},n,lr_multiplier,0.05,m,f_update,updaterule,optimizeall);
        end
        
        %% fc -> peak_filter_gradient_scale = 0.2*fs/nfft, shelving_filter_gradient_scale = 0.1*fs/nfft
        if m == 1 || m == numlayers
            [fc{m},mfc{m},vfc{m},gfc{m}] = optimize(fc{m},dydfc{m},mfc{m},vfc{m},gfc{m},n,lr_multiplier,0.1*fs/1024,m,f_update,updaterule,optimizeall);%0.5
        else
            [fc{m},mfc{m},vfc{m},gfc{m}] = optimize(fc{m},dydfc{m},mfc{m},vfc{m},gfc{m},n,lr_multiplier,0.2*fs/1024,m,f_update,updaterule,optimizeall);%2
        end
        
        %% fb -> peak_filter_gradient_scale = 2*fs/nfft
        if m~=1 && m~= numlayers
            [fb{m},mfb{m},vfb{m},gfb{m}] = optimize(fb{m},dydfb{m},mfb{m},vfb{m},gfb{m},n,lr_multiplier,2*fs/1024,m,f_update,updaterule,optimizeall);%1
        end
    end
end

% copy utility %
copy2struct;

end

    