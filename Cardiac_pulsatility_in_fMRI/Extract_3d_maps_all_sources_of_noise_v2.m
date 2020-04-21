% Correct the name of this file !!

clear, clc, close all

baseDir_PC1='../../RawData/';
baseDir='\\DiskStation\HCP\';

load([baseDir_PC1,'Subject_list.mat'])
subject_list = [subject_list_R3];
subj_remove = {'199251','114924'};
for i = 1:length(subj_remove)
    ind = find(subject_list == subj_remove{i} );  subject_list(ind) = [];
end
subject_list = subject_list(1:100);
nSubj = length(subject_list); nScans = nSubj*4;
task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};


save_path = ['Export/2020_02_01/Stat_maps_CP_crossValid/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

TR = 0.72;
HPF_f = 0.008;    [filt_b,filt_a] = butter(2,HPF_f*2*TR,'high');

nii_MNI = load_untouch_nii('E:\CloudStation\Atlases\standard\MNI152_T1_2mm_brain_mask.nii.gz'); imgMask = nii_MNI.img;
nii_1 = load_untouch_nii('E:\CloudStation\Atlases\MNI152_2mm_1vol.nii.gz');
nii_3 = load_untouch_nii('E:\CloudStation\Atlases\MNI152_2mm_3vol.nii.gz');
nii_12 = load_untouch_nii('E:\CloudStation\Atlases\MNI152_2mm_12vol.nii.gz');
% nii_6 = load_untouch_nii('E:\CloudStation\Atlases\MNI152_2mm_6vol.nii.gz');

volDel=20; volDel_from_end=0; NV_0 = 1200;

kFold=3; NV = NV_0-volDel-volDel_from_end;
block_length=round((NV)/kFold); ind_blocks=zeros(kFold,2);
for i=1:kFold
    ind_blocks(i,1)=1+(i-1)*block_length;
    ind_blocks(i,2)=min(i*block_length,NV);  ind_blocks(kFold,2) = NV;
end

%% GLM - Motion, Respiration, Cardiac (RETROICOR), Stimulus


for c = 1 : nScans
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    path_output=[save_path,'/',subject,'_',task];
    
    baseDir_subj = [baseDir,'RawData_R3_131/',subject,'/',task];
    filepath_nii=[baseDir_subj,'/rfMRI_',task,'.nii.gz'];
    %         filepath_nii=[baseDir_subj,'/func.nii.gz'];
    
    filepath_mask=[baseDir_subj,'/brainmask_fs.2.nii.gz'];
    filepath_movRegr=[baseDir_subj,'/Movement_Regressors_dt.txt'];
        
    filepath_MRacq=[baseDir_PC1,'Physio/',subject,'_',task,'/phys.mat']; 
%     load(filepath_MRacq,'Fs','trig','TR'); Ts = 1/Fs;
    load(filepath_MRacq); Ts = 1/Fs;
    load([baseDir_PC1,'Physio/',subject,'_',task,'/Phys_sum.mat']); 
    
            %% ----------------------------------------     
        disp('Loading data'), tic
        
        nii=load_untouch_nii(filepath_nii); img=single(nii.img);
        nii_mask=load_untouch_nii(filepath_mask); imgMask=nii_mask.img;
        TR = 0.72;
        
        volDel=20;                img(:,:,:,1:volDel) = [];
        [NX,NY,NZ,NV]=size(img);
        
        ind_BOLD=find(trig==1);     trig(ind_BOLD(1:volDel))=0;    ind_BOLD=find(trig==1);
        time = 0:Ts:(length(trig)-1)*Ts;        
        timeMR=time(find(trig==1));     
        
        disp('Loading data - done'), fprintf('Time elapsed : %3.1f  minutes \n\n ', toc/60), 
    
    
        %% ================================
    disp('Extracting  img_col '), tic
    
    NVX=0;
    CordAll=zeros(NX*NY*NZ,3);
    for x=1:NX
        for y=1:NY
            for z=1:NZ
                if imgMask(x,y,z)>0.5
                    voxel=reshape(img(x,y,z,:),NV,1);
                    voxelMean = mean(voxel);
                    if voxelMean>100
                        NVX=NVX+1;
                        CordAll(NVX,1)=x;
                        CordAll(NVX,2)=y;
                        CordAll(NVX,3)=z;
                        voxel=single(voxel); voxel=detrend(voxel,'linear');
                        img(x,y,z,:)=voxel+voxelMean;
                    end
                end
            end
        end
    end
    CordAll=CordAll(1:NVX,:);
    
    img_col=zeros(NV,NVX);
    for vox=1:NVX
        x=CordAll(vox,1);    y=CordAll(vox,2);    z=CordAll(vox,3);
        img_col(:,vox)=reshape(img(x,y,z,:),NV,1);
    end
    img_col = filtfilt(filt_b,filt_a,img_col);
    
    disp('Extracting  img_col  - done'), fprintf('Time elapsed: %3.1f minutes \n\n ', toc/60)
    
    figure
    GS = mean(img_col,2);
    plot( timeMR, GS), title('Global signal (GS)'), xlim([min(timeMR),max(timeMR)])
    
    
    %% Framewise displacement    
    
    movRegr=load(filepath_movRegr);        movRegr=[movRegr, movRegr.^2];
    Mov = abs(movRegr(:,7:12));    Mov1 = sum(Mov(:,1:3),2) ;    Mov2 =  50 *  sum(Mov(:,4:6),2) * (2*pi()/360) ;    FD = Mov1 + Mov2;
    FD(1:volDel) = []; movRegr(1:volDel,:) = [];    
    movRegr = filtfilt(filt_b,filt_a,movRegr);
    
    fprintf('Mean FD: %3.2f mm \n', mean(FD))
    
    figure
    plot(timeMR, zscore(FD))   
    
    %%  Disentagle resp, movRegr and DVARS
    
%     resp_der = func_RETR_Resp(...)    % Load regressors for RETROICOR    (resp)
    
    regr_nuis = [regr_resp, movRegr, GS];    
    regr_nuis = filtfilt(filt_b,filt_a,regr_nuis ) ;
    regr_nuis = [regr_nuis, ones(NV,1)];
    
    parfor vox=1:NVX  
        voxel = img_col(:,vox);                
        B = regr_nuis\voxel;    
        img_col(:,vox) = voxel - regr_nuis*B;        
    end
        
    
    %% Cardiac pulsatility - RETROICOR
    
    M_order = 6 ;
    RETR_CardRegr_cont=RETR_Card_regressors_v2(time,PPGlocs,M_order);
    RETR_CardRegr = RETR_CardRegr_cont(ind_BOLD+round(+0.4*Fs),:);
    RETR_CardRegr = filtfilt(filt_b,filt_a,RETR_CardRegr);
    for j = 1:size(RETR_CardRegr,2)
        voxel = RETR_CardRegr(:,j);
        B = regr_nuis\voxel;
        RETR_CardRegr(:,j) = voxel - regr_nuis*B;
    end
    regr_RETR = [RETR_CardRegr, ones(NV,1)];  
    
    %% Cardiac pulsatility - DCPM

    HR =HRV;
    M_order = 6; N = length(time);
    memory = 60/mean(HR);
    CPM_IR = func_CPM_cos(Ts, memory, M_order);
    
    u = zeros(size(time)); uA = u;
    nPeaks = length(PPGlocs);
    for i = 1:nPeaks
        t = PPGlocs(i);
        [val loc] = min(abs(time-t));
        u(loc) = 1;      uA(loc) = cardiac(loc);
    end
    
    CPM_regr_all = zeros(N+length(CPM_IR)-1,M_order);
    CPM_Amp_regr_all = zeros(N+length(CPM_IR)-1,M_order);
    for m = 1:M_order*2
        x = conv(u,CPM_IR(:,m));    CPM_regr_all(:,m) = x;
        x = conv(uA,CPM_IR(:,m));  CPM_Amp_regr_all(:,m) = x;
    end
    
    shift = +0.9;
    CPM_regr_all = CPM_regr_all(ind_BOLD+round(shift*Fs),:);
    CPM_regr_all = filtfilt(filt_b,filt_a,CPM_regr_all);    
    for j = 1:size(CPM_regr_all,2)
        voxel = CPM_regr_all(:,j);
        B = regr_nuis\voxel;
        CPM_regr_all(:,j) = voxel - regr_nuis*B;
    end
    CPM_regr_all = [CPM_regr_all, ones(NV,1)];
    
    CPM_Amp_regr_all = CPM_Amp_regr_all(ind_BOLD+round(shift*Fs),:);
    CPM_Amp_regr_all = filtfilt(filt_b,filt_a,CPM_Amp_regr_all);    
    for j = 1:size(CPM_Amp_regr_all,2)
        voxel = CPM_Amp_regr_all(:,j);
        B = regr_nuis\voxel;
        CPM_Amp_regr_all(:,j) = voxel - regr_nuis*B;
    end
    CPM_Amp_regr_all = [CPM_Amp_regr_all, ones(NV,1)];    
    
    
    
    %%  Extract 3d stat maps and export them   -----------------------------
    
    disp('Extracting 3d stat maps - RETROICOR '), tic
    close all
            
    R_all_col = zeros(NVX,3);
    R_noCV_col = zeros(NVX,1);
    B_CPM_CA_col = zeros(NVX,12);
    parfor vox=1:NVX  % NVX
        voxel = img_col(:,vox);
        
        r_valid_k=zeros(3,kFold);
            for k=1:kFold
                ind_valid = ind_blocks(k,1):ind_blocks(k,2);
                ind_train = 1:NV; ind_train(ind_valid)=[];
                voxel_train =  voxel(ind_train); voxel_valid = voxel(ind_valid);
                
                regr = regr_RETR;
                regr_train = regr(ind_train,:);  regr_valid = regr(ind_valid,:);
                B = regr_train\voxel_train;     yPred = regr_valid*B;
                r_valid_k(1,k) = corr(voxel_valid,yPred)   ;
                
                regr = CPM_regr_all;
                regr_train = regr(ind_train,:);  regr_valid = regr(ind_valid,:);
                B = regr_train\voxel_train;     yPred = regr_valid*B;
                r_valid_k(2,k) = corr(voxel_valid,yPred)   ;

                regr = CPM_Amp_regr_all;
                regr_train = regr(ind_train,:);  regr_valid = regr(ind_valid,:);
                B = regr_train\voxel_train;     yPred = regr_valid*B;
                r_valid_k(3,k) = corr(voxel_valid,yPred)   ;
            end              
        R_all_col(vox,:) = mean(r_valid_k,2) ;
        
        B = CPM_regr_all\voxel;     yPred = CPM_regr_all*B;
        R_noCV_col(vox,:) = corr(voxel,yPred)  ;
        B_CPM_CA_col(vox,:) = B(1:12);                
    end
    
    mean(R_all_col,1)
    
    R_all = zeros(NX,NY,NZ,3);
    R_noCV_all = zeros(NX,NY,NZ);
    B_CPM_CA = zeros(NX,NY,NZ,12);
    for vox=1:NVX
        x=CordAll(vox,1);    y=CordAll(vox,2);    z=CordAll(vox,3);
        R_all(x,y,z,:) = R_all_col(vox,:);
        R_noCV_all(x,y,z) = R_noCV_col(vox,:);
        B_CPM_CA(x,y,z,:) = B_CPM_CA_col(vox,:);
    end
        
    nii_3.img=R_all;         save_untouch_nii(nii_3,[path_output,'_CP_R_all.nii.gz'])  
    nii_1.img=R_noCV_all;         save_untouch_nii(nii_1,[path_output,'_CP_R_CPM_CA_noCV.nii.gz'])  
    nii_12.img=B_CPM_CA;         save_untouch_nii(nii_12,[path_output,'_CP_B_CPM_CA.nii.gz'])  
          
    %  --------------------------------------    
    
    close all
    disp('Extracting  3d stat maps  - done'), fprintf('Time elapsed: %3.1f  \n ', toc/60),        
    
    
    
    %%   Extract B parameters
    
    
    
    
    
end


















