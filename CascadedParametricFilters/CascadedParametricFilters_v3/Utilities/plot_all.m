%% impulse Response
% figure(3);
% subplot(2,1,1);
% stem(yl_impulse); grid on;
% xlabel('samples','Interpreter','latex'); ylabel('amplitude','Interpreter','latex');
% title('Desired IR','Interpreter','latex');
% subplot(2,1,2);
% stem(yl_hat_impulse(1:n)); grid on;
% xlabel('samples','Interpreter','latex'); ylabel('amplitude','Interpreter','latex');
% title('Approximation IR','Interpreter','latex');
% subplot(2,1,1);
% plot(20*log10(abs(yl_impulse)+eps)); grid on; hold on;
% xlabel('samples','Interpreter','latex'); ylabel('amplitude','Interpreter','latex');
% title('Desired IR','Interpreter','latex'); 
% subplot(2,1,2);
% plot(20*log10(abs(yl_hat_impulse(1:n))+eps)); grid on;
% xlabel('samples','Interpreter','latex'); ylabel('amplitude','Interpreter','latex');
% title('Approximation IR','Interpreter','latex');hold off;

%% magnitude response (IR)
YL                  = 20*log10(abs(fft(yl_impulse,NFFT)));
YL                  = YL(1:NFFT/2+1);
YL_HAT              = 20*log10(abs(fft(yl_hat_impulse,NFFT)));
YL_HAT              = YL_HAT(1:NFFT/2+1);
YL_INIT             = 20*log10(abs(fft(yl_hat_impulse_init,NFFT)));
YL_INIT             = YL_INIT(1:NFFT/2+1);
YL_DESIRED_SMOOTH   = H_smooth;

if plotresults
    %% approximation, target, initial with fc lines
    figure(4); 
    semilogx(fvec,YL_HAT,'LineWidth',3,'Color',hsured); hold on
    semilogx(fvec,YL_DESIRED_SMOOTH,'LineWidth',3,'Color',hsublack);
    semilogx(fvec,YL_INIT,'LineWidth',2,'Color',hsucyan);
    grid on;
    xlabel('Frequency in Hz','Interpreter','latex');ylabel('Magnitude in dB','Interpreter','latex');
    if H_smooth(1)<0
        legend('Approximation','Target','Initial','Interpreter','latex','AutoUpdate', 'off','location','northwest');
    else
        legend('Approximation','Target','Initial','Interpreter','latex','AutoUpdate', 'off','location','southwest');
    end
    xlim([20 20000]);
    for i = 1:numel(paramstructbest.fc)-1
        plot([paramstructbest.fc{i} paramstructbest.fc{i}],ylim,'Color',rgbval,'LineWidth',0.5);
        plot([paramstructinit.fc{i} paramstructinit.fc{i}],ylim,'Color',rgbval*0,'LineWidth',0.5);
    end
    hold off;
    set(gca,'fontname','times');
    set(findall(gcf,'-property','FontSize'),'FontSize',40);
    set(gca,'linewidth',2.5);
    
    %% approximation, target, initial, error without fc lines
    figure(5);
    semilogx(fvec,YL_HAT,'LineWidth',3,'Color',hsured);hold on
    semilogx(fvec,YL_DESIRED_SMOOTH,'LineWidth',3,'Color',hsublack);
    semilogx(fvec,YL_INIT,'LineWidth',2,'Color',hsucyan);
    semilogx(fvec,YL_DESIRED_SMOOTH-YL_HAT,'LineWidth',2,'Color',hsuredlight);
    grid on;
    xlabel('Frequency in Hz','Interpreter','latex');ylabel('Magnitude in dB','Interpreter','latex');
    if H_smooth(1)<0
        legend('Approximation','Target','Initial','Error','Interpreter','latex','AutoUpdate', 'off','location','northwest');
    else
        legend('Approximation','Target','Initial','Error','Interpreter','latex','AutoUpdate', 'off','location','southwest');
    end
    xlim([20 20000]);
    hold off;
    set(gca,'fontname','times');
    set(findall(gcf,'-property','FontSize'),'FontSize',40);
    set(gca,'linewidth',2.5);
    
    %% approximation, original non-smooth target without fc lines
    figure(6);
    semilogx(fvec,YL_HAT,'LineWidth',3,'Color',hsured); hold on
    semilogx(fvec,YL,'LineWidth',3,'Color',hsublack);
    grid on;
    xlabel('Frequency in Hz','Interpreter','latex');ylabel('Magnitude in dB','Interpreter','latex');
    if H_smooth(1)<0
        legend('Approximation','Target','Interpreter','latex','AutoUpdate', 'off','location','northwest');
    else
        legend('Approximation','Target','Interpreter','latex','AutoUpdate', 'off','location','southwest');
    end
    xlim([20 20000]);
    hold off;
    set(gca,'fontname','times');
    set(findall(gcf,'-property','FontSize'),'FontSize',40);
    set(gca,'linewidth',2.5);
    
    %% approximation, target, error without fc lines
    figure(7);
    semilogx(fvec,YL_HAT,'LineWidth',3,'Color',hsured);hold on
    semilogx(fvec,YL_DESIRED_SMOOTH,'LineWidth',3,'Color',hsublack);
    semilogx(fvec,YL_DESIRED_SMOOTH-YL_HAT,'LineWidth',2,'Color',hsuredlight);
    grid on;
    xlabel('Frequency in Hz','Interpreter','latex');ylabel('Magnitude in dB','Interpreter','latex');
    if H_smooth(1)<0
        legend('Approximation','Target','Error','Interpreter','latex','AutoUpdate', 'off','location','northwest');
    else
        legend('Approximation','Target','Error','Interpreter','latex','AutoUpdate', 'off','location','southwest');
    end
    xlim([20 20000]);
    hold off;
    set(gca,'fontname','times');
    set(findall(gcf,'-property','FontSize'),'FontSize',40);
    set(gca,'linewidth',2.5);
    
    %% error vs epoch
    figure(8);
    plot(error,'b');
    ylabel('MLSD in dB','Interpreter','latex'); xlabel('Iterations','Interpreter','latex'); grid on;
    title('error vs iterations','Interpreter','latex');
    set(gca,'fontname','times');
    set(findall(gcf,'-property','FontSize'),'FontSize',40);
    set(gca,'linewidth',2.5);
    
    %% approximation vs initial, filters
    figure(9);
    [s1,s4,m1,m2] = fcn_visualize(paramstructinit,paramstructbest,length(yl_impulse),fs,[]);
    grid on;
    xlabel('Frequency in Hz','Interpreter','latex');ylabel('Magnitude in dB','Interpreter','latex');
    if H_smooth(1)<0
        legend([s1,s4],'Approximation','Initial','Interpreter','latex','AutoUpdate', 'off','location','northwest');
    else
        legend([s1,s4],'Approximation','Initial','Interpreter','latex','AutoUpdate', 'off','location','southwest');
    end
    xlim([50 20000]);
    for i = 1:numel(paramstructbest.fc)-1
        plot([paramstructbest.fc{i} paramstructbest.fc{i}],ylim,'Color',hsured,'LineWidth',1);
        plot([paramstructinit.fc{i} paramstructinit.fc{i}],ylim,'Color',hsucyan,'LineWidth',1);
    end
    set(gca,'fontname','times');
    set(findall(gcf,'-property','FontSize'),'FontSize',40);
    set(gca,'linewidth',2.5);
    hold off;
    
    %% calculate metrices
    E_init      = mean(abs(YL_DESIRED_SMOOTH(2:464)-YL_INIT(2:464)));
    E_init_max  = max(abs(YL_DESIRED_SMOOTH(2:464)-YL_INIT(2:464)));
    E_init_min  = min(abs(YL_DESIRED_SMOOTH(2:464)-YL_INIT(2:464)));
    E_fin       = mean(abs(YL_DESIRED_SMOOTH(2:464)-YL_HAT(2:464)));
    E_fin_max   = max(abs(YL_DESIRED_SMOOTH(2:464)-YL_HAT(2:464)));
    E_fin_min   = min(abs(YL_DESIRED_SMOOTH(2:464)-YL_HAT(2:464)));
    N_filt      = numel(paramstructbest.fc);
    max(E_fin_max)
    
    %% save
%     filename = ['Folder_to_save_in/',subj_name,'_',num2str(az_angle),'_',num2str(el_angle),'.mat'];
%     save(filename,'yl_impulse','yl_hat_impulse','yl_hat_impulse_init','H_smooth','error','paramstructbest','paramstructinit','fs');
else
    %% calculate metrices 
    E_init = mean(abs(YL_DESIRED_SMOOTH(2:464)-YL_INIT(2:464)));
    E_init_max = max(abs(YL_DESIRED_SMOOTH(2:464)-YL_INIT(2:464)));
    E_init_min = min(abs(YL_DESIRED_SMOOTH(2:464)-YL_INIT(2:464)));
    E_fin = mean(abs(YL_DESIRED_SMOOTH(2:464)-YL_HAT(2:464)));
    E_fin_max = max(abs(YL_DESIRED_SMOOTH(2:464)-YL_HAT(2:464)));
    E_fin_min = min(abs(YL_DESIRED_SMOOTH(2:464)-YL_HAT(2:464)));
    N_filt = numel(paramstructbest.fc);
end