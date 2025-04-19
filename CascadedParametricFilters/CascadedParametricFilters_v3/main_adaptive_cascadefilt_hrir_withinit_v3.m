close all
clear
clc

addpath('Initialization','Modules','Utilities','../10Peaks');

%% color templates

rgbval      = [0.5 0.5 0.5];
hsured      = [0.77255 0 0.25882];
hsugray     = [0.5647 0.5176 0.4627];
hsublack    = [0 0 0];
hsuredlight = [216 150 170]./256;
hsucyan     = [0,0.4431,0.4157];

%% only plot without saving

plotresults = 1;

%% network dependent variables

updaterule  = 'adam';
numepochs   = 20;
optimizeall = 1;

%% variables related to filter initialization

N_smooth        = 12;
onlypeaks       = 0;
ITD             = 0;
NFFT            = 1024;

%% dataset selection

fs              = 44100;
az_angle_vec    = [0, 55, 80];    % -55
el_angle_vec    = [-45 0 45 90];
load_path       = '../10Peaks/subject_';
subject_list    = {'008'};   
n_ear           = 1;        % 1: Left Ear, 2: Right Ear

for sub = 1 : size(subject_list)
    for az = 1 : numel(az_angle_vec)
        for el = 1 : numel(el_angle_vec)
                
            subj_name   = subject_list{sub};    
            az_angle  	= az_angle_vec(az);
            el_angle  	= el_angle_vec(el);

            % load data
            load([load_path,subj_name,'/approximation_az',num2str(az_angle),'_el',num2str(el_angle),'.mat'],...
                'h_original','param_approx','H_mean');
            h_original      = h_original(n_ear,:);
            param_approx    = param_approx(:,:,n_ear);
            mean_coarse     = H_mean(n_ear);
            clear H_mean
            
            % extend IR, smooth, initialize
            preprocess_approximated_hrtfs;
            
            % adapt network parameters
            lr_multiplier   = 8e-2;
            for k = 1 : numepochs
                tic
                for m = 1 : numfilters % option to optimize 1 filter at a time or all together
                        e_plot                          = [];
                        paramstruct.y{1}                = x; % input
                        [y,paramstruct,global_G,e_plot] = adaptive_cascadefilt(yl, yi, fs, paramstruct, numfilters, lr_multiplier,global_G,e_plot,m,updaterule,NFFT,H_smooth_label,optimizeall,k,numepochs);
                        yi                              = y{end};
                        
                        %% Store the error values for monitoring %%%%%%%%%%
                        error       = [error e_plot];
                        errorepoch  = [errorepoch mean(e_plot)];
                        
                        %% Save the parametricmodel with the lowest error %
                        if k == 1
                            errorbest       = mean(e_plot);
                            paramstructbest = paramstruct;
                        end
                        if k > 1
                            if errorepoch(end) < errorepoch(end-1)
                                errorbest = errorepoch(end);
                            else
                                lr_multiplier = lr_multiplier*0.96;
                            end
                            if errorepoch(end) <= errorbest
                                paramstructbest = paramstruct;
                            end
                        end
                end
                if n_ear == 1
                    save(['Optimization/subject_',subj_name,'/left_ear_az',num2str(az_angle),'_el',num2str(el_angle),'.mat'],...
                        'paramstructinit','paramstructbest','H_coarse','H_smooth','mean_coarse','error','k','numepochs','fs');
                else
                    save(['Optimization/subject_',subj_name,'/right_ear_az',num2str(az_angle),'_el',num2str(el_angle),'.mat'],...
                        'paramstructinit','paramstructbest','H_coarse','H_smooth','mean_coarse','error','k','numepochs','fs');
                end
                toc
                ['Epoch #',num2str(k),' of ',num2str(numepochs),' done']
            end
            
            % impulse response
            Length                  = K * n;
            xn                      = [1 zeros(1,Length-1)];
            paramstructlabel.y{1}   = xn;
            paramstructinit.y{1}    = xn;
            yl_impulse              = s;
            paramstructbest.y{1}    = xn;
            yl_hat_impulse          = adaptive_cascadefiltforward(fs, paramstructbest, numfilters);
            yl_hat_impulse_init     = adaptive_cascadefiltforward(fs, paramstructinit, numfilters);
            
            % plot results
            plot_all;
            % close all
        end
    end
end

