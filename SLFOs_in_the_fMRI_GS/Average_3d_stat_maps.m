

clear, clc, close all

baseDir ='../../RawData/';
load_path = ['Export/2020_03_24/PPG_correlation_maps/'];
path_output = ['Export/2020_03_24/PPG_correlation_maps_group/'];  if exist(path_output,'file')~=7, mkdir(path_output); end


load([baseDir,'Subject_list.mat'])

subject_list = [subject_list_R3];
subj_remove = {'199251','114924'};
for i = 1:length(subj_remove)
    ind = find(subject_list == subj_remove{i} );  subject_list(ind) = [];
end
subject_list = subject_list(1:100);
n_subj = length(subject_list); nScans = n_subj*4;
task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};

nii_MNI = load_untouch_nii('E:\CloudStation\Atlases\MNI152_2mm_1vol.nii.gz');
img_MNI =single(nii_MNI.img);
[ NX, NY, NZ ]   =  size(img_MNI);

%% GLM - PPG, PPGconv

r_PPG_all = zeros(n_subj, NX,NY,NZ);
r_PPGconv_all = zeros(n_subj, NX,NY,NZ);

tic
parfor s =   1 :  n_subj
    subject = char(subject_list(s));
    fprintf('Subject: %s     (%d/%d)  \n',subject,s,n_subj)    
            
    for task_ID = 1
        task = char(task_list(task_ID));
        fprintf('  ------------------------------------------------------- \n ')
        fprintf('Task: %s      (%d/%d)  \n',task,task_ID,4)
        
        path_output_subj = [load_path,subject,'_',task]; 
        
        nii_subj = load_untouch_nii([path_output_subj,'_r_PPG.nii.gz']);
        r_img_subj = single(nii_subj.img);
        r_PPG_all(s,:,:,:) = r_img_subj;
        
        nii_subj = load_untouch_nii([path_output_subj,'_r_PPGconv.nii.gz']);
        r_img_subj = single(nii_subj.img);
        r_PPGconv_all(s,:,:,:) = r_img_subj;                
    end    
end

fprintf('Time elapsed (minutes): %3.1f  \n', toc/60),

%%   ----------------------------------------

nii_MNI.img = squeeze(mean(r_PPG_all,1)); save_untouch_nii(nii_MNI,[path_output,'r_PPG.nii.gz'])
nii_MNI.img = squeeze(mean(r_PPGconv_all,1)); save_untouch_nii(nii_MNI,[path_output,'r_PPGconv.nii.gz'])

load chirp,  sound(y,Fs)





    



