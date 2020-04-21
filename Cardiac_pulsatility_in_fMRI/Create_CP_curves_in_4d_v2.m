

clear, clc, close all

baseDir ='../../RawData/';


load_path = ['Export/2020_02_01/Stat_maps_CP_crossValid/'];
path_output = ['Export/2020_02_01/CP_curves_statistics/'];  if exist(path_output,'file')~=7, mkdir(path_output); end

load([baseDir,'Subject_list.mat'])

subject_list = [subject_list_R3];
subj_remove = {'199251','114924'};
for i = 1:length(subj_remove)
    ind = find(subject_list == subj_remove{i} );  subject_list(ind) = [];
end
subject_list = subject_list(1:100);
n_subj = length(subject_list); nScans = n_subj*4;
task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};


nii_stats_4vol=load_untouch_nii('Export/2019_11_01/MNI152_2mm_4vol.nii.gz' );
nii_stats_41vol=load_untouch_nii('Export/2019_11_01/MNI152_2mm_41vol.nii.gz' );
nii_MNI = load_untouch_nii('Export/2019_11_01/MNI152_2mm_1vol.nii.gz');
img_MNI =single(nii_MNI.img);
[ NX, NY, NZ ]   =  size(img_MNI);

%% GLM - Motion, Respiration, Cardiac (RETROICOR), Stimulus

t_IR = linspace(0,1,41);      CP_cycle = 1;  N_IR = length(t_IR); Ts = t_IR(2)-t_IR(1);

% RETR_basis =  func_CPM_cos(Ts, CP_cycle, 4);
DCPM_basis =  func_CPM_cos(Ts, CP_cycle, 6);


r_RETR_LR_CV_all = zeros(n_subj, NX,NY,NZ);
r_DCPM_LR_CV_all = zeros(n_subj, NX,NY,NZ);
r_RETR_RL_CV_all = zeros(n_subj, NX,NY,NZ);
r_DCPM_RL_CV_all = zeros(n_subj, NX,NY,NZ);

r_DCPM_LR_all = zeros(n_subj, NX,NY,NZ);
r_DCPM_RL_all = zeros(n_subj, NX,NY,NZ);

CP_curves_DCPM_LR_all = zeros(NX,NY,NZ, N_IR);
CP_curves_DCPM_RL_all = zeros(NX,NY,NZ, N_IR);

tic
parfor s =   1 :  n_subj
    subject = char(subject_list(s));
    fprintf('Subject: %s     (%d/%d)  \n',subject,s,n_subj)    
    
        
    for task_ID = 1
        task = char(task_list(task_ID));
        fprintf('  ------------------------------------------------------- \n ')
        fprintf('Task: %s      (%d/%d)  \n',task,task_ID,4)
        
        path_output_subj = [load_path,subject,'_',task]; 
        
        nii_subj = load_untouch_nii([path_output_subj,'_CP_R_all.nii.gz']);
        r_img_subj = single(nii_subj.img);
        r_RETR_LR_CV_all(s,:,:,:) = squeeze(r_img_subj(:,:,:,1));
        r_DCPM_LR_CV_all(s,:,:,:) = squeeze(r_img_subj(:,:,:,2));        
        
        nii_subj = load_untouch_nii([path_output_subj,'_CP_R_CPM_CA_noCV.nii.gz']);
        r_img_subj = single(nii_subj.img);
        r_DCPM_LR_all(s,:,:,:) = r_img_subj;        
        
        nii_subj = load_untouch_nii([path_output_subj,'_CP_B_CPM_CA.nii.gz']);
        img_subj = single(nii_subj.img);
        ind_vox = find(r_img_subj>0.01);
        img4d = zeros(NX,NY,NZ,N_IR);
        for v = 1 : length(ind_vox)
            vox = ind_vox(v);            [x, y, z] = ind2sub(size(img_subj), vox);
            B_vox = squeeze(img_subj(x, y, z, :));
            IR =  B_vox' * DCPM_basis'; IR = IR(:)./max(abs(IR));      %        plot(IR),    r_img_subj(x,y,z)  , hold on
            img4d(x,y,z,:) = IR  *  r_img_subj(x,y,z)   ;
        end
        CP_curves_DCPM_LR_all = CP_curves_DCPM_LR_all +  img4d/n_subj  ;
    end
    
    for task_ID = 2
        task = char(task_list(task_ID));
        fprintf('  ------------------------------------------------------- \n ')
        fprintf('Task: %s      (%d/%d)  \n',task,task_ID,4)
        
        path_output_subj = [load_path,subject,'_',task];        
                
        nii_subj = load_untouch_nii([path_output_subj,'_CP_R_all.nii.gz']);
        r_img_subj = single(nii_subj.img);
        r_RETR_RL_CV_all(s,:,:,:) = squeeze(r_img_subj(:,:,:,1));
        r_DCPM_RL_CV_all(s,:,:,:) = squeeze(r_img_subj(:,:,:,2));        
        
        nii_subj = load_untouch_nii([path_output_subj,'_CP_R_CPM_CA_noCV.nii.gz']);
        r_img_subj = single(nii_subj.img);
        r_DCPM_RL_all(s,:,:,:) = r_img_subj;        
        
        nii_subj = load_untouch_nii([path_output_subj,'_CP_B_CPM_CA.nii.gz']);
        img_subj = single(nii_subj.img);
        ind_vox = find(r_img_subj>0.01);
        img4d = zeros(NX,NY,NZ,N_IR);
        for v = 1 : length(ind_vox)
            vox = ind_vox(v);            [x, y, z] = ind2sub(size(img_subj), vox);
            B_vox = squeeze(img_subj(x, y, z, :));
            IR =  B_vox' * DCPM_basis'; IR = IR(:)./max(abs(IR));      %        plot(IR),    r_img_subj(x,y,z)  , hold on
            img4d(x,y,z,:) = IR  *  r_img_subj(x,y,z)   ;
        end
        CP_curves_DCPM_RL_all = CP_curves_DCPM_RL_all +  img4d/n_subj  ;
    end
    
end

fprintf('Time elapsed (minutes): %3.1f  \n', toc/60),

%%   ----------------------------------------

r_RETR_LR_CV_all = zeros(n_subj, NX,NY,NZ);
r_DCPM_LR_CV_all = zeros(n_subj, NX,NY,NZ);
r_RETR_RL_CV_all = zeros(n_subj, NX,NY,NZ);
r_DCPM_RL_CV_all = zeros(n_subj, NX,NY,NZ);

r_DCPM_LR_all = zeros(n_subj, NX,NY,NZ);
r_DCPM_RL_all = zeros(n_subj, NX,NY,NZ);


r_RETR_LR_CV_mean = squeeze(mean(r_RETR_LR_CV_all,1));
nii_MNI.img = r_RETR_LR_CV_mean; save_untouch_nii(nii_MNI,[path_output,'r_RETR_LR_CV_mean.nii.gz'])

r_DCPM_LR_CV_mean = squeeze(mean(r_DCPM_LR_CV_all,1));
nii_MNI.img = r_DCPM_LR_CV_mean; save_untouch_nii(nii_MNI,[path_output,'r_DCPM_LR_CV_all.nii.gz'])

r_RETR_RL_CV_mean = squeeze(mean(r_RETR_RL_CV_all,1));
nii_MNI.img = r_RETR_RL_CV_mean; save_untouch_nii(nii_MNI,[path_output,'r_RETR_RL_CV_mean.nii.gz'])

r_DCPM_RL_CV_mean = squeeze(mean(r_DCPM_RL_CV_all,1));
nii_MNI.img = r_DCPM_RL_CV_mean; save_untouch_nii(nii_MNI,[path_output,'r_DCPM_RL_CV_all.nii.gz'])

r_DCPM_RL_CV_mean = squeeze(mean(r_DCPM_LR_all,1));
nii_MNI.img = r_DCPM_RL_CV_mean; save_untouch_nii(nii_MNI,[path_output,'r_DCPM_LR_noCV_all.nii.gz'])

r_DCPM_RL_CV_mean = squeeze(mean(r_DCPM_RL_all,1));
nii_MNI.img = r_DCPM_RL_CV_mean; save_untouch_nii(nii_MNI,[path_output,'r_DCPM_RL_noCV_all.nii.gz'])

%  T-test across models   ------------------------

ttest_RETR_to_DCPM = zeros(NX,NY,NZ);
ttest_DCPM_LR_to_RL = zeros(NX,NY,NZ);
parfor x = 1:NX
    for y = 1:NY
        for z = 1:NZ
            dist_1 =  squeeze(r_RETR_LR_CV_all(:,x,y,z));
            dist_2 =  squeeze(r_DCPM_LR_CV_all(:,x,y,z));
            [h,p,ci,stats] = ttest(dist_2,dist_1);
            t = stats.tstat;            t(isnan(t))=0;
            ttest_RETR_to_DCPM(x,y,z) = t;
            
            dist_1 =  squeeze(r_DCPM_LR_all(:,x,y,z));
            dist_2 =  squeeze(r_DCPM_RL_all(:,x,y,z));
            [h,p,ci,stats] = ttest(dist_2,dist_1);
            t = stats.tstat;            t(isnan(t))=0;
            ttest_DCPM_LR_to_RL(x,y,z) = t;            
        end
    end
end

nii_MNI.img = ttest_RETR_to_DCPM; save_untouch_nii(nii_MNI,[path_output,'ttest_RETR_to_DCPM_CV.nii.gz'])
nii_MNI.img = ttest_DCPM_LR_to_RL; save_untouch_nii(nii_MNI,[path_output,'ttest_DCPM_LR_to_RL_noCV.nii.gz'])


%  Curves   ------------------------
nii_stats_41vol.img = CP_curves_DCPM_LR_all; save_untouch_nii(nii_stats_41vol,[path_output,'CPM_CA_curves_LR.nii.gz'])
nii_stats_41vol.img = CP_curves_DCPM_RL_all; save_untouch_nii(nii_stats_41vol,[path_output,'CPM_CA_curves_RL.nii.gz'])
CP_curves_LR_m_RL = CP_curves_DCPM_RL_all - CP_curves_DCPM_LR_all;
nii_stats_41vol.img = CP_curves_LR_m_RL; save_untouch_nii(nii_stats_41vol,[path_output,'CPM_CA_curves_LR_to_RL.nii.gz'])



load chirp,  sound(y,Fs)


%%  ----------------------------------------
%  Save maps for a specific subject

if 0
    CP_curves_DCPM_all = CP_curves_DCPM_all * n_subj/4;
    % r_LBMM_mean = squeeze(r_LBMM_all(:);
    % nii_MNI.img = r_LBMM_mean; save_untouch_nii(nii_MNI,[path_output,'sbj05_R_CPM_CA_LR.nii.gz'])
    nii_stats_41vol.img = CP_curves_DCPM_all; save_untouch_nii(nii_stats_41vol,[path_output,'sbj79_CPM_CA_curves_LR.nii.gz'])
    load chirp,  sound(y,Fs)
end
%%   ----------------------------------------



tmp_L = zeros(41,1);
tmp_R = zeros(41,1);

for x = 1:NX
    for y = 1:NY
        for z = 1:NZ
            tmp = squeeze(CP_curves_DCPM_all(x,y,z,:));
            if y <round(NY/2)
                tmp_L = tmp_L + tmp;
            else
                tmp_R = tmp_R + tmp;                
            end
        end
    end
end

figure
plot(tmp_L), hold on
plot(tmp_R)



flag_kMeans = 0;
if flag_kMeans
    %%  --------------------------
    % keep only voxels prone toa rtifacts ...
    
    r_col = reshape(r_LBMM_mean,[NX*NY*NZ,1]);
    
    [ ind ]  = find(r_col>0.25);
    tmp = reshape(CP_curves_DCPM_all, [NX*NY*NZ,41]);
    tmp = tmp(ind,:);
    
    k = 4;
    [idx,C,sumd,D] = kmeans(tmp,k,'Distance','correlation');
    
    C = C';
    Cnorm = zeros(size(C));
    for i =1:k
        x = C(:,i);
        x = x - x(1); x = x/ max(abs(x));
        Cnorm(:,i) = x
    end    
    plot(Cnorm)
        
    %%
    
    id_col = zeros(NX*NY*NZ,1);
    id_col(ind) = idx;
    
    id_3d = reshape(id_col, [NX,NY,NZ]);
    nii_MNI.img = id_3d; save_untouch_nii(nii_MNI,[path_output,'id_3d_thr_0p25.nii.gz'])
    
    
    
end




