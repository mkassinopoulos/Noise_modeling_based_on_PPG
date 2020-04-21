
function [Ts,Fs,TR, img_GM_Card_col,trig,PPGlocs, HR, BR, RF, resp, DVARS, cardiac, GS, movRegr] = load_scan(subject,task,baseDir,filepath_MRacq,flag_order)

if flag_order==1
%     load(['Export/2019_10_15/ROIs_w_cardiac_artifacts/',subject,'_',task]);
    load(['Export/2020_01_27/ROIs_w_cardiac_artifacts/',subject,'_',task]);
    img_GM_Card_col = img_WB_Card_col;
elseif flag_order==2
    load(['Export/2019_10_15/ROIs_w_breathing_artifacts/',subject,'_',task]);
    img_GM_Card_col = img_GM_Breath_col;
elseif flag_order == 0
%         load(['Export/2019_10_15/ROIs_w_randomVoxels/',subject,'_',task]);
            load(['Export/2020_01_27/ROIs_w_cardiac_artifacts/',subject,'_',task]);
                img_GM_Card_col = img_WB_rand_col;
elseif flag_order == 3
        load(['Export/2019_10_15/ROIs_w_random_voxels_breathing/',subject,'_',task]);
            img_GM_Card_col = img_GM_Breath_col;
elseif flag_order == 5
    img_GM_Card_col = [];
end
load([baseDir,'/Physio/',subject,'_',task,'/Phys_sum.mat']);
load(filepath_MRacq,'Fs','trig','TR','cardiac'); Ts = 1/Fs;

filepath_movRegr=[baseDir,'/Physio/',subject,'_',task,'/Movement_Regressors_dt.txt'];
movRegr=load(filepath_movRegr); 


BR = BR(:);
HR = HRV(:);
RF = IHF(:);

filepath_input=[baseDir,'Atlas/',subject,'_',task,'/'];
load([filepath_input,'TissueBasedRegressors_1199.mat'])

GS = WB.MA;

end




