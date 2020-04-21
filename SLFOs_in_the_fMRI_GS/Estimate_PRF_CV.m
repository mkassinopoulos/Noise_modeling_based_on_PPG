
clear, close all, clc


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

save_path = ['Export/2020_03_24/Estimate_PRF_shift_RV_CV/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

shift = 5;     % Lag time for PPG-Amp

%%  1:  Load physiological variables (heart rate and respiration) and global signal (GS) from MAT-File

%  Set the following parameters !!

Ts = 1/400;
% HPF_f = 0.7;  [filt_b,filt_a] = butter(2,HPF_f*2*Ts,'high');

TR = 0.72;
HPF_f = 0.008;  [filt_b,filt_a] = butter(2,HPF_f*2*TR,'high');

volDel=30; NV_0 = 1200; NV = NV_0 - volDel;
Ts_10 = 0.1; Fs_10 = 10; t_win= 0 :Ts_10:60; N_win = length(t_win);

kFold=3; block_length=round((NV)/kFold); ind_blocks=zeros(kFold,2);
for i=1:kFold
    ind_blocks(i,1)=1+(i-1)*block_length;
    ind_blocks(i,2)=min(i*block_length,NV);  ind_blocks(kFold,2) = NV;
end

r_all = zeros(nScans,3);  tic
r_CV_stand = zeros(nScans,1); r_CV_PPG = zeros(nScans,1);
r_CV_PPG_only_PPG = zeros(nScans,1);

parfor c = 1 :    nScans
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    [Ts,Fs,TR, trig,PPGlocs, HR, resp_10, DVARS, cardiac, GS, FD, ~, ~, movRegr]  =  load_scan(subject,task,baseDir);
    
    HRmean(c)  = mean(HR);     HRstd(c) = std(HR);    FDmean(c) = mean(FD);
    GSmean(c) = mean(GS);     GSstd(c) = std(GS);
    
    time = 0:Ts:(length(trig)-1)*Ts;
    time_10 = 0:0.1:time(end);
    ind_BOLD=find(trig==1);
    trig(ind_BOLD(1:volDel))=0;
    ind_BOLD=find(trig==1);
    timeMR=time(trig==1); timeMR = timeMR +TR/2;
    NV = length(timeMR);
    
    GS(1:volDel) = [];      GS = filtfilt(filt_b,filt_a,GS ) ;
    GS = zscore(GS);
    
    ind_BOLD_10 = zeros(NV,1);
    for i = 1:NV
        t = timeMR(i);        [val loc] = min(abs(time_10-t));
        ind_BOLD_10(i) = loc;
    end
    
    %% Extract RV
    
    RV = zeros(size(resp_10));
    N = length(resp_10);
    for i = 2:N-1
        ind_1 = i-6*Fs_10;   ind_1 = max(1,ind_1);
        ind_2 = i+6*Fs_10;   ind_2 = min(ind_2, N);
        RV(i) = std(resp_10(ind_1:ind_2));
    end
    RV(1) = RV(2); RV(end) = RV(end-1);    RV = RV(:);
    
    %% 2: Estimate PRF parameters
    
    resp_10 = filloutliers(resp_10, 'clip','ThresholdFactor',4);    
    resp_s = smooth(resp_10,10*1.5) ;
    RF=diff(resp_s); RF=[0;RF(:)]; RF = RF.^2;
    
    RF = RV;    
%     HR = cardiac;
    
    ga_opts = gaoptimset('TolFun',1e-10,'StallGenLimit',20,'Generations',100,'Display','none','UseParallel',1);   % Display: iter
    options = optimoptions('fmincon','Display','off','Algorithm','interior-point',...
        'UseParallel',true,'MaxIterations',100,'MaxFunctionEvaluations',3000,'OptimalityTolerance',1e-8);    % 'PlotFcn','optimplotfval'
    
    PRF_par = [  3.1    2.5   5.6    0.9    1.9   2.9   12.5    0.5    3.1    2.5   5.6    0.9  ]; PRF_par_0 = PRF_par;
    ub = PRF_par+3;
    lb = PRF_par-3; lb(find(lb<0))=0;
    
    
    % -------------------------------
    
    r_valid_k_stand=zeros(kFold,1);
    r_valid_k_PPG=zeros(kFold,1);
    r_valid_k_PPG_onlyPPG=zeros(kFold,1);
    for k=1:kFold
        ind_valid = ind_blocks(k,1):ind_blocks(k,2);
        ind_train = 1:NV; ind_train(ind_valid)=[];
        if kFold==1, ind_train=1:NV; end
        
        PRF_par = PRF_par_0;
        % Standard  -----------------------------
        h_train = @(P) func_PRF_sc(P,Ts_10,HR,RF,ind_BOLD_10,GS,1, ind_train, ind_valid,filt_b,filt_a);
        %                 PRF_par = ga(h_train,length(ub),[],[],[],[],lb,ub,[],[],ga_opts);
%         PRF_par = fmincon(h_train,PRF_par,[],[],[],[],lb,ub,[],options);
        
        h = @(P) func_PRF_sc(P,Ts_10,HR,RF,ind_BOLD_10,GS, 0, ind_train, ind_valid,filt_b,filt_a);
        [obj_function,CRF_sc,RRF_sc,HR_conv,RF_conv,r_PRF_sc,yPred, HR_conv_MR, RF_conv_MR] = h(PRF_par);
        r_valid_k_stand(k) = r_PRF_sc(1);
        
        PRF_par = PRF_par_0;
        % with PPG  -----------------------------
%         h_train = @(P) func_PRF_w_ePPG_NoConv(P,Ts_10,HR,RF,cardiac, ind_BOLD_10,GS,1,ind_train,ind_valid,filt_b,filt_a);
        h_train = @(P) func_PRF_w_ePPG_shift(P,Ts_10,HR,RF,cardiac, ind_BOLD_10,GS,1,ind_train,ind_valid,filt_b,filt_a,shift);
%         h_train = @(P) func_PRF_w_ePPG(P,Ts_10,HR,RF,cardiac, ind_BOLD_10,GS,1,ind_train,ind_valid,filt_b,filt_a);
%         h_train = @(P) func_PRF_w_ePPG_noConv(P,Ts_10,HR,RF,cardiac, ind_BOLD_10,GS,1,ind_train,ind_valid,filt_b,filt_a);
        %                 PRF_par = ga(h_train,length(ub),[],[],[],[],lb,ub,[],[],ga_opts);
        PRF_par = fmincon(h_train,PRF_par,[],[],[],[],lb,ub,[],options);
        
%         h = @(P) func_PRF_w_ePPG(P,Ts_10,HR,RF,cardiac,ind_BOLD_10,GS,0,ind_train,ind_valid,filt_b,filt_a);
%         h = @(P) func_PRF_w_ePPG_noConv(P,Ts_10,HR,RF,cardiac,ind_BOLD_10,GS,0,ind_train,ind_valid,filt_b,filt_a);
%         [obj_function,CRF_sc,RRF_sc,HR_conv,RF_conv,PPG_conv, r_PRF_sc,yPred, HR_conv_MR, RF_conv_MR] = h(PRF_par);

        h = @(P) func_PRF_w_ePPG_shift(P,Ts_10,HR,RF,cardiac,ind_BOLD_10,GS,0,ind_train,ind_valid,filt_b,filt_a,shift);
        [obj_function,CRF_sc,RRF_sc,BRF_sc, r_PRF_sc] = h(PRF_par);
        
%         h = @(P) func_PRF_w_ePPG_NoConv(P,Ts_10,HR,RF,cardiac,ind_BOLD_10,GS,0,ind_train,ind_valid,filt_b,filt_a);
%         [obj_function,CRF_sc,RRF_sc, r_PRF_sc] = h(PRF_par);

        r_valid_k_PPG(k) = r_PRF_sc(1);
        r_valid_k_PPG_onlyPPG(k) = r_PRF_sc(4);
        
    end
    r_CV_stand(c) = mean(r_valid_k_stand);
    r_CV_PPG(c) = mean(r_valid_k_PPG);
    r_CV_PPG_only_PPG(c) = mean(r_valid_k_PPG_onlyPPG);
    
    %     h_train = @(P) func_PRF_w_ePPG(P,Ts_10,HR,RF,ind_BOLD_10,GS,1,uePPG);
    %     % PRF_par = ga(h_train,length(ub),[],[],[],[],lb,ub,[],[],ga_opts);
    %     PRF_par = fmincon(h_train,PRF_par,[],[],[],[],lb,ub,[],options);
    %
    %     h = @(P) func_PRF_w_ePPG(P,Ts_10,HR,RF,ind_BOLD_10,GS,0,uePPG);
    %     [obj_function,CRF_sc,RRF_sc,HR_conv,RF_conv,r_PRF_sc,yPred, HR_conv_MR, RF_conv_MR] = h(PRF_par);
    
    
    fprintf(' ----------------------------------------------- \n')
    fprintf('Correlation b/w GS and PRF output \n')
    %     fprintf('CRF (HR): %3.2f  \n',r_PRF_sc(2))
    %     fprintf('RRF (RF): %3.2f  \n',r_PRF_sc(3))
    %     fprintf('CRF & RRF (HR & RF): %3.2f  \n',r_PRF_sc(1))
    fprintf('HR+RF -- cross-validation: %3.2f  \n',mean(r_valid_k_stand))
    fprintf('HR+RF+PPG -- cross-validation: %3.2f  \n',mean(r_valid_k_PPG))
    fprintf('HR+RF+PPG -- only PPGconv: %3.2f  \n',mean(r_valid_k_PPG_onlyPPG))
    
    
    
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60), load chirp,  sound(y,Fs)

save([save_path, 'RV_PPG_shift_5_ext'],'r_CV_stand','r_CV_PPG','r_CV_PPG_only_PPG')

%%  ---------------------------------

plot(r_all)

fprintf('CV - stand:  %3.2f  \n', mean(r_CV_stand))
fprintf('CV - RF+HR+ePPG:  %3.2f  \n', mean(r_CV_PPG))

plot(r_CV_PPG-r_CV_stand)

[ttest_h ttest_p] = ttest(r_CV_PPG,r_CV_stand)


x1_sbj = reshape(r_CV_stand',[4 100]);  x1_sbj = mean(x1_sbj)';
x1_sbj = reshape(r_CV_PPG',[4 100]);  x1_sbj = mean(x1_sbj)';
x2_sbj = reshape(r_CV_PPG',[4 100]);  x2_sbj = mean(x2_sbj)';

fprintf('subject-based: CV - stand:  %3.2f+-%3.2f  \n', mean(x1_sbj),std(x1_sbj))
fprintf('subject-based: CV - RF+HR+ePPG:  %3.2f+-%3.2f  \n', mean(x2_sbj),std(x2_sbj))

[ttest_p ttest_h] = ttest(x1_sbj,x2_sbj)


%%    Compare corr(GS,PA) vs corr(GS,PAconv)
% Load workspace..

if 0   

mean(r_GS_PPG_NoConv)
mean(r_GS_PPG_conv)

x1_sbj = reshape(r_GS_PPG_NoConv',[4 100]);  x1_sbj = abs(mean(x1_sbj)');
x2_sbj = reshape(r_GS_PPG_conv',[4 100]);  x2_sbj = mean(x2_sbj)';

fprintf('subject-based: CV - stand:  %3.2f+-%3.2f  \n', mean(x1_sbj),std(x1_sbj))
fprintf('subject-based: CV - RF+HR+ePPG:  %3.2f+-%3.2f  \n', mean(x2_sbj),std(x2_sbj))

[ttest_p ttest_h] = ttest(x1_sbj,x2_sbj)

end



