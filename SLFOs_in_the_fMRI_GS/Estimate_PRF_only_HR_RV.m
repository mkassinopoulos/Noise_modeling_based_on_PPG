
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


save_path = ['Export/2020_02_01/Estimate_PRF_only_HR_RV/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

%%  1:  Load physiological variables (heart rate and respiration) and global signal (GS) from MAT-File

%  Set the following parameters !!

Ts = 1/400; 
% HPF_f = 0.7;  [filt_b,filt_a] = butter(2,HPF_f*2*Ts,'high');

TR = 0.72; 
HPF_f = 0.01;  [filt_b,filt_a] = butter(2,HPF_f*2*TR,'high');

volDel=30;
Ts_10 = 0.1; Fs_10 = 10; t_win= 0 :Ts_10:60; N_win = length(t_win);


r_all = zeros(nScans,3);  tic
CRF_all = zeros(N_win,nScans);
RRF_all = zeros(N_win,nScans);
FC_all = zeros(nScans, 3,3);
FC_partial_all = zeros(nScans,3,3);
%%  -----------

for c = 1   : nScans
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
    
    % Extract RV
    
    RV = zeros(size(resp_10));
    N = length(resp_10);    
    for i = 2:N-1
        ind_1 = i-6*Fs_10;   ind_1 = max(1,ind_1);
        ind_2 = i+6*Fs_10;   ind_2 = min(ind_2, N);
        RV(i) = std(resp_10(ind_1:ind_2));
    end
    RV(1) = RV(2); RV(end) = RV(end-1);    RV = RV(:);
    
    %% 2: Estimate PRF parameters
            
    RF = RV;
    %     HR = HR./uePPG;
    
    ga_opts = gaoptimset('TolFun',1e-10,'StallGenLimit',20,'Generations',100,'Display','iter','UseParallel',1);   % Display: iter
    options = optimoptions('fmincon','Display','off','Algorithm','interior-point',...
        'UseParallel',true,'MaxIterations',100,'MaxFunctionEvaluations',3000,'OptimalityTolerance',1e-8);    % 'PlotFcn','optimplotfval'
    
    PRF_par = [  3.1    2.5   5.6    0.9    1.9   2.9   12.5    0.5    ]; PRF_par_0 = PRF_par;
    ub = PRF_par+3;
    lb = PRF_par-3; lb(find(lb<0))=0;
    
    
        PRF_par = PRF_par_0;
        % Standard  -----------------------------
        h_train = @(P) func_PRF_sc(P,Ts_10,HR,RF,ind_BOLD_10,GS,1, 1:NV, 1:NV,filt_b,filt_a);
        %                 PRF_par = ga(h_train,length(ub),[],[],[],[],lb,ub,[],[],ga_opts);
        PRF_par = fmincon(h_train,PRF_par,[],[],[],[],lb,ub,[],options);
        
        h = @(P) func_PRF_sc(P,Ts_10,HR,RF,ind_BOLD_10,GS, 0, 1:NV, 1:NV,filt_b,filt_a);
        [obj_function,CRF_sc,RRF_sc,HR_conv,RF_conv,r_PRF_sc,yPred, HR_conv_MR, RF_conv_MR] = h(PRF_par);
        
    
    CRF_all(:,c) = CRF_sc;
    RRF_all(:,c) = RRF_sc;       
    r_all(c,:) = r_PRF_sc;

%     X = [yPred_resp, yPred_card,GS];
%     FC = corr(X),    
%     FC_partial = partialcorr(X);        
%     FC_all(c,:,:) = FC;        
%     FC_partial_all(c,:,:) = FC_partial;
    
    
    fprintf(' ----------------------------------------------- \n')
    fprintf('Correlation b/w GS and PRF output \n')
    fprintf('CRF (HR): %3.2f  \n',r_PRF_sc(2))
    fprintf('RRF (RF): %3.2f  \n',r_PRF_sc(3))
    fprintf('CRF & RRF (HR & RF): %3.2f  \n',r_PRF_sc(1))
    %  plot(BRF_sc)
         
    
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60), load chirp,  sound(y,Fs)

save([save_path, 'only_HR_RV'],'r_all','CRF_all','RRF_all')

%%  ---------------------------------

plot(r_all)

mean(r_all)

CM = squeeze(mean(FC_all))
imagesc(CM)


CM_partial = squeeze(mean(FC_partial_all))
imagesc(CM_partial)


mean(r_all_wos)
mean(r_all_ws)

x1 = r_all_wos(:,1);
x2 = r_all_ws(:,1);

x1_sbj = reshape(x1',[4 100]);  x1_sbj = mean(x1_sbj)';
x2_sbj = reshape(x2',[4 100]);  x2_sbj = mean(x2_sbj)';

[ttest_p ttest_h] = ttest(x1_sbj,x2_sbj)


%%  --------------------

nScans = 400;

CRF_mean = zeros(N_win,1);
BRF_mean = zeros(N_win,1);
RRF_mean = zeros(N_win,1);
for c = 1:nScans
    x = CRF_all(:,c)*r_all(c,2);
    CRF_mean = CRF_mean + x/nScans;
    
    x = BRF_all(:,c)*r_all(c,4);
    BRF_mean = BRF_mean + x/nScans;
    
    x = RRF_all(:,c)*r_all(c,3);
    RRF_mean = RRF_mean + x/nScans;
end

figure('position', [  940   703   555   420])
ax = plot(t_win,[CRF_mean, BRF_mean, RRF_mean],'linewidth',3), grid on

ax = gca;
ax.GridLineStyle = '--';
ax.Box = 'off';
ax.XAxisLocation = 'origin';

legend('CRF','PARF','RRF');  legend boxoff

xlim([0 50]), ylim([-0.5 0.5])
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')







