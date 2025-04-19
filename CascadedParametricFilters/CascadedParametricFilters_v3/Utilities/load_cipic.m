function [h_left,h_right,ITD] = load_cipic(cipic_path, subj_name, az_angle, el_angle)

    load(strcat(cipic_path,subj_name,'/hrir_final.mat'),'hrir_l','hrir_r','ITD');

    az_cipic    = [-80, -65, -55, -45:5:45, 55, 65, 80];
    el_cipic    = -45:5.625:230.625;

    az_indx     = find(az_cipic == az_angle);
    el_indx     = find(el_cipic == el_angle);

    h_left      = squeeze(hrir_l(az_indx,el_indx,:));
    h_right     = squeeze(hrir_r(az_indx,el_indx,:));
    
    ITD         = ITD(az_indx,el_indx);
    
end