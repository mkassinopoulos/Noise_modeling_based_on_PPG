

clear, clc, close all

baseDir='../../RawData/';
load([baseDir,'Subject_list.mat'])
subject_list = [subject_list_R3];
subj_remove = {'199251','114924'};
for i = 1:length(subj_remove)
    ind = find(subject_list == subj_remove{i} );  subject_list(ind) = [];
end
subject_list = subject_list(1:100);
nSubj = length(subject_list); nScans = nSubj*4;
task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};


M_order = 8;    % 8th order RETROICOR  --> 16 DoF
nSV = 1000;

volDel=20; volDel_from_end=10;

kFold=3; NV = 1180-volDel_from_end;
block_length=round((NV)/kFold); ind_blocks=zeros(kFold,2);
for i=1:kFold
    ind_blocks(i,1)=1+(i-1)*block_length;
    ind_blocks(i,2)=min(i*block_length,NV);  ind_blocks(kFold,2) = NV;
end
temp_shift_scale = -2:0.1:2; nShifts = length(temp_shift_scale);

TR = 0.72;
HPF_f = 0.008;    [filt_b,filt_a] = butter(2,HPF_f*2*TR,'high');


flag_order = 0 ;
for model_card = 1:3    % 1: RETROICOR, 2: DCPM_CA, 3: DCPM_VA

if flag_order ==1
    txt_1 = '_Cardiac'
else
    txt_1 = '_random'
end

if model_card == 1
    txt_2 = '_RETROICOR'
elseif model_card == 2
    txt_2 = '_DCPM_CA'
elseif model_card == 3
    txt_2 = '_DCPM_VA'
end

save_path = ['Export/2020_02_01/CP_Model_order_selection/'];  if exist(save_path,'file')~=7, mkdir(save_path); end
save_filename = ['r_all',txt_1,txt_2];


%%  --------------------------------

tic
r_all = zeros(nScans,nShifts, nSV,M_order);
HRmean_all = zeros(nScans,1);
HRstd_all = zeros(nScans,1);
parfor c = 1 : nScans
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    filepath_MRacq=[baseDir,'/Physio/',subject,'_',task,'/phys.mat'];
    [Ts,Fs,TR, img_GM_Card_col,trig,PPGlocs, HR, ~, ~, resp, ~, cardiac, GS, movRegr] =  load_scan(subject,task,baseDir,filepath_MRacq, flag_order);
    
    HRmean = mean(HR);
    HRstd = std(HR);
    HRmean_all(c) = HRmean;
    HRstd_all(c) = HRstd;
    
    data = img_GM_Card_col(1:end-volDel_from_end,:);
    data =  filtfilt(filt_b,filt_a,data);
    N = length(cardiac);

      %%  ---------------------------------------------
    
    ind_BOLD=find(trig==1);  
    trig(ind_BOLD(1:volDel))=0;    trig(ind_BOLD(end-volDel_from_end+1:end)) = 0;    
    ind_BOLD=find(trig==1);
    time = 0:Ts:(length(trig)-1)*Ts;
    timeMR=time(trig==1);
    NV = length(timeMR);
       
    
    %%  Regressors for BM, HM and GS  ---------------------------------
    
    GS(1:volDel) = []; GS = GS(1:end-volDel_from_end);
    ind = volDel+1: volDel+NV;
    movRegr = movRegr(ind,1:12);
    
    resp_der = diff(resp); resp_der = [resp_der(1); resp_der];
    regr_resp = [resp, resp_der, resp.^2, resp_der.^2];
    regr_resp = regr_resp(ind_BOLD,:);
    
    regr_nuis = [regr_resp, movRegr, GS];
    regr_nuis = zscore(filtfilt(filt_b,filt_a,regr_nuis )) ;
    regr_nuis = [regr_nuis, [1:NV]', ones(NV,1)];
    
    for v = 1:nSV
        voxel = data(:,v);
        B = regr_nuis\voxel;
        data(:,v) = voxel - regr_nuis*B;
    end
    
    %%  ------------------------
    
    memory = 60/mean(HR);
    CPM_IR = func_CPM_cos(Ts, memory, M_order);
    
    u = zeros(size(time)); uA = u;
    nPeaks = length(PPGlocs);
    for i = 1:nPeaks
        t = PPGlocs(i);
        [val loc] = min(abs(time-t));
        u(loc) = 1;      uA(loc) = cardiac(loc);
    end
    
    CPM_regr_all = zeros(N,M_order);
    CPM_Amp_regr_all = zeros(N,M_order);
    for m = 1:M_order*2
        x = conv(u,CPM_IR(:,m));    CPM_regr_all(:,m) = x(1:N);
        x = conv(uA,CPM_IR(:,m));  CPM_Amp_regr_all(:,m) = x(1:N);
    end
    
    RETR_regr_all = RETR_Card_regressors_v2(time,PPGlocs,M_order);
    
    switch model_card
        case 1
            regr_all_x = RETR_regr_all;
        case 2
            regr_all_x = CPM_regr_all;
        case 3
            regr_all_x = CPM_Amp_regr_all;
    end
        
    %%  -------------------------------
        
    r_temp_shift = zeros(nShifts,nSV,M_order);
    for  c_temp = 1:nShifts
        shift = temp_shift_scale(c_temp);
        RETR_CardRegr = regr_all_x(ind_BOLD+round(shift*Fs),:);
        RETR_CardRegr = filtfilt(filt_b,filt_a,RETR_CardRegr);
        
        for j = 1:size(RETR_CardRegr,2)
            voxel = RETR_CardRegr(:,j);
            B = regr_nuis\voxel;
            RETR_CardRegr(:,j) = voxel - regr_nuis*B;
        end        
        regr = [RETR_CardRegr, ones(NV,1)];
        
        %%  ---------------------------------------------
        
        r_scan = zeros(nSV,M_order);
        for v = 1:nSV
            voxel = data(:,v);            
            r_valid_k=zeros(M_order,kFold);
            for k=1:kFold
                ind_valid = ind_blocks(k,1):ind_blocks(k,2);
                ind_train = 1:NV; ind_train(ind_valid)=[];
                if kFold==1, ind_train=1:NV; ind_valid = 1:NV; end
                
                regr_train = regr(ind_train,:);  regr_valid = regr(ind_valid,:);
                voxel_train =  voxel(ind_train); voxel_valid = voxel(ind_valid);                                
                for m = 1:M_order
                    ind_regr = [1:2*m,size(regr,2)];
                    B = regr_train(:,ind_regr)\voxel_train;     yPred = regr_valid(:,ind_regr)*B;
                    r_valid_k(m,k) = corr(voxel_valid,yPred)   ;            
                end
                
            end
            r_scan(v,:)=mean(r_valid_k,2);
        end
        r_temp_shift(c_temp,:,:) = r_scan;
    end
    
%     tmp = squeeze(mean(r_temp_shift,2));    max(tmp(:))
%     plot(temp_shift_scale,tmp), grid on
    
%     tmp_RETR = tmp(:,4);
%     tmp_DCPM_VA = tmp(:,4);
%     plot(temp_shift_scale, tmp_RETR), hold on
%     plot(temp_shift_scale, tmp_DCPM_VA), hold on
    
    r_all(c,:,:,:) = r_temp_shift;
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60),
load chirp,  sound(y,Fs)


save([save_path,save_filename],'r_all','temp_shift_scale');   % -v7.3

end






