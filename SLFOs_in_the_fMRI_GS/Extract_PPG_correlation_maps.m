

clear, clc, close all

% baseDir='../../RawData/';
% addpath(baseDir)
% load('subject_list.mat');
% subject_list = [subject_list_R1; subject_list_R2];

baseDir='\\DiskStation\HCP\';

% save_path = ['Export/2019_10_15/ROIs_w_randomVoxels/'];  if exist(save_path,'file')~=7, mkdir(save_path); end
save_path = ['Export/2020_03_24/PPG_correlation_maps/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};
load([baseDir,'RawData_PC1_Michalis/Subject_list.mat'])

subject_list = [subject_list_R3];
subj_remove = {'199251','114924'};
for i = 1:length(subj_remove)
    ind = find(subject_list == subj_remove{i} );  subject_list(ind) = [];
end
subject_list = subject_list(1:100);
n_subj = length(subject_list); nScans = n_subj*4;
task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};

HPF_f = 0.008; TR = 0.72;
[filt_b,filt_a] = butter(2,HPF_f*2*TR,'high');

nii_MNI = load_untouch_nii('E:\CloudStation\Atlases\MNI152_2mm_1vol.nii.gz');
img_MNI = nii_MNI.img;


%% GLM - Motion, Respiration, Cardiac (RETROICOR), Stimulus


for s = 1 : n_subj
    subject = char(subject_list(s));
    fprintf('  ===================================== \n\n ')
    fprintf('Subject: %s     (%d/%d)  \n\n',subject,s,n_subj)
    %     path_sub = [path_init,'/',subject_ID,'/anat/'];
    
    for task_ID = 1   % :4
        close all
        
        task = char(task_list(task_ID));
        fprintf('  ------------------------------------------------------- \n ')
        fprintf('Task: %s      (%d/%d)  \n\n',task,task_ID,4)
        
        path_output=[save_path,'/',subject,'_',task];
        
        baseDir_subj = [baseDir,'RawData_R3_131/',subject,'/',task];
        filepath_nii=[baseDir_subj,'/rfMRI_',task,'.nii.gz'];
        %         filepath_nii=[baseDir_subj,'/func.nii.gz'];
        
        filepath_mask=[baseDir_subj,'/brainmask_fs.2.nii.gz'];
        filepath_movRegr=[baseDir_subj,'/Movement_Regressors_dt.txt'];
        
        filepath_MRacq=[baseDir,'RawData_PC1_Michalis/Physio/',subject,'_',task,'/phys.mat'];
        %         load(filepath_MRacq,'Fs','trig','TR'); Ts = 1/Fs;
        load(filepath_MRacq); Ts = 1/Fs;
        load([baseDir,'RawData_PC1_Michalis/Physio/',subject,'_',task,'/Phys_sum.mat']);
        
        
        %% ----------------------------------------
        
        disp('Loading data'), tic
        nii=load_untouch_nii(filepath_nii); img=single(nii.img);
        nii_mask=load_untouch_nii(filepath_mask); imgMask=nii_mask.img;
        
        TR = 0.72;
        
        volDel=30;                img(:,:,:,1:volDel) = [];
        [NX,NY,NZ,NV]=size(img);
        
        ind_BOLD=find(trig==1);     trig(ind_BOLD(1:volDel))=0;    ind_BOLD=find(trig==1);
        time = 0:Ts:(length(trig)-1)*Ts;
        
        timeMR=time(find(trig==1));
        
        ind_BOLD_10 = zeros(NV,1);
        for i = 1:NV
            t = timeMR(i);        [val loc] = min(abs(time_10-t));
            ind_BOLD_10(i) = loc;
        end
        
        disp('Loading data - done'), fprintf('Time elapsed : %3.1f  minutes \n\n ', toc/60),
        
        %% ================================
        disp('Extracting  img_col and tissue maps '), tic
        
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
        
        img_col=zeros(NV,NVX);     Tissue_col = zeros(NVX,1); thresh_cardiac_col = zeros(NVX,1);
        for vox=1:NVX
            x=CordAll(vox,1);    y=CordAll(vox,2);    z=CordAll(vox,3);
            img_col(:,vox)=reshape(img(x,y,z,:),NV,1);
        end
        ind_cardiac = find(thresh_cardiac_col==1);
        img_col = filtfilt(filt_b,filt_a,img_col);
        
        disp('Extracting  img_col and tissue maps - done'), fprintf('Time elapsed: %3.1f minutes \n\n ', toc/60)
        
        GS =mean(img_col,2);
        
        figure
        plot( zscore(smooth(GS,10))), title('Global signal (GS; smoothed)'), xlim([min(timeMR),max(timeMR)]), hold on
        
        
        %%   Load RETROICOR regressors
        
        RETR_RespRegr = [];
        load([baseDir,'RawData_PC1_Michalis/Physio/',subject,'_',task,'/RETROICOR.mat']);
        
        RETR_RespRegr = RETR_RespRegr(volDel+1:end,1:4);
        
        %% Cardiac pulsatility - DCPM
        
        HR =HRV;
        M_order = 6; N = length(time);
        memory = 60/mean(HR);
        CPM_IR = func_CPM_cos(Ts, memory, M_order);
        
        uA = zeros(size(time));
        nPeaks = length(PPGlocs);
        for i = 1:nPeaks
            t = PPGlocs(i);
            [val loc] = min(abs(time-t));
            uA(loc) = cardiac(loc);
        end
        
        CPM_Amp_regr_all = zeros(N+length(CPM_IR)-1,M_order);
        for m = 1:M_order*2
            x = conv(uA,CPM_IR(:,m));  CPM_Amp_regr_all(:,m) = x;
        end
        shift = +0.9;       CPM_Amp_regr_all = CPM_Amp_regr_all(ind_BOLD+round(shift*Fs),:);
        
        movRegr=load(filepath_movRegr);         movRegr(1:volDel,:) = [];
        
        regr = [movRegr, RETR_RespRegr, CPM_Amp_regr_all];
        regr = filtfilt(filt_b,filt_a,regr); regr = detrend(regr,'linear');
        regr = [regr, ones(NV,1)];
        
        
        %%  Extract PPG-Amp
        
        indPPGlocs = zeros(size(PPGlocs));
        for i = 1 :length(PPGlocs)
            [val loc] = min(abs(PPGlocs(i)-time));
            indPPGlocs(i) = loc;
        end
        pks = cardiac(indPPGlocs);
        cardiac = interp1([0,time(indPPGlocs),time(end)],[pks(1),pks',pks(end)],time);
        cardiac = interp1(time, cardiac, time_10);    cardiac = zscore(cardiac(:));
        
        
        %% Extract RV
        
        Fs_10  = 10; Ts_10 = 0.1;
        RV = zeros(size(resp_10));
        N = length(resp_10);
        for i = 2:N-1
            ind_1 = i-6*Fs_10;   ind_1 = max(1,ind_1);
            ind_2 = i+6*Fs_10;   ind_2 = min(ind_2, N);
            RV(i) = std(resp_10(ind_1:ind_2));
        end
        RV(1) = RV(2); RV(end) = RV(end-1);    RV = RV(:);
        
        %% 2: Estimate PRF parameters
        
        
        ga_opts = gaoptimset('TolFun',1e-10,'StallGenLimit',20,'Generations',100,'Display','iter','UseParallel',1);   % Display: iter
        options = optimoptions('fmincon','Display','off','Algorithm','interior-point',...
            'UseParallel',true,'MaxIterations',100,'MaxFunctionEvaluations',3000,'OptimalityTolerance',1e-8);    % 'PlotFcn','optimplotfval'
        
        PRF_par = [  3.1    2.5   5.6    0.9    1.9   2.9   12.5    0.5    3.1    2.5   5.6    0.9  ]; PRF_par_0 = PRF_par;
        ub = PRF_par+3;
        lb = PRF_par-3; lb(find(lb<0))=0;
        
        h_train = @(P) func_PRF_w_ePPG_tmp(P,Ts_10,HR,RV,cardiac, ind_BOLD_10,GS,1,1:NV,1:NV,filt_b,filt_a);
        %                 PRF_par = ga(h_train,length(ub),[],[],[],[],lb,ub,[],[],ga_opts);
        PRF_par = fmincon(h_train,PRF_par,[],[],[],[],lb,ub,[],options);
        
        h = @(P) func_PRF_w_ePPG_tmp(P,Ts_10,HR,RV,cardiac,ind_BOLD_10,GS,0,1:NV,1:NV,filt_b,filt_a);
        %     [obj_function,CRF_sc,RRF_sc,BRF_sc,HR_conv,RF_conv,r_PRF_sc,yPred, HR_conv_MR, RF_conv_MR, yPred_PPG] = h(PRF_par);
        
        [obj_function,CRF_sc,RRF_sc,BRF_sc, r_PRF_sc,yPred, yPred_card, yPred_resp, yPred_PPG] = h(PRF_par);
        
        
        fprintf(' ----------------------------------------------- \n')
        fprintf('Correlation b/w GS and PRF output \n')
        fprintf('CRF (HR): %3.2f  \n',r_PRF_sc(2))
        fprintf('RRF (RF): %3.2f  \n',r_PRF_sc(3))
        fprintf('PPG (BP): %3.2f  \n',r_PRF_sc(4))
        fprintf('CRF & RRF (HR & RF): %3.2f  \n',r_PRF_sc(1))
        
        
        %%  ----------------------------
        
        r_PPG_col = zeros(NVX,1);
        r_PPGconv_col = zeros(NV,1);
        parfor i = 1:NVX
            voxel = img_col(:,i);                B=regr\voxel;   yPred=regr*B;
            voxel_clean = voxel-yPred;
            
            r_PPG_col(i) = corr(voxel_clean, cardiac(ind_BOLD_10));
            r_PPGconv_col(i) = corr(voxel_clean, yPred_PPG);
        end
        
        %%  ---------------------
        
                
        r_PPG = zeros(NX,NY,NZ);
        r_PPGconv = zeros(NX,NY,NZ);
        for vox=1:NVX
            x=CordAll(vox,1);    y=CordAll(vox,2);    z=CordAll(vox,3);
            r_PPG(x,y,z,:) = r_PPG_col(vox,:);
            r_PPGconv(x,y,z,:) = r_PPGconv_col(vox,:);
        end
        
        nii_MNI.img=r_PPG;         save_untouch_nii(nii_MNI,[path_output,'_r_PPG.nii.gz'])
        nii_MNI.img=r_PPGconv;         save_untouch_nii(nii_MNI,[path_output,'_r_PPGconv.nii.gz'])
        
        
        close all
        
    end
end


















