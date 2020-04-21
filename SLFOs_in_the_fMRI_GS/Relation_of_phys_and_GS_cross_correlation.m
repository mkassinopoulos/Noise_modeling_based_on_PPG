
%%    Clear workspace

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

save_path = ['Export/2020_02_01/Relation_of_phys_and_GS_CrossCorrelation/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

path_filename = 'workspace_phys_interactions';

%  --------------------------------

Fs_10 = 10;  Fs = 400; Ts_10 = 0.1;
r_all = zeros(nScans,2);

HRmean = zeros(nScans,1);
HRstd = zeros(nScans,1);
FDmean = zeros(nScans,1);
GSmean = zeros(nScans,1);
GSstd = zeros(nScans,1);

N_IR = 1001;

Br_to_HR_IR = zeros(nScans,N_IR);
Br_to_PPG_IR = zeros(nScans,N_IR);
Br_to_GS_IR = zeros(nScans,N_IR);
HR_to_Br_IR = zeros(nScans,N_IR);
HR_to_PPG_IR = zeros(nScans,N_IR);
HR_to_GS_IR = zeros(nScans,N_IR);
PPG_to_Br_IR = zeros(nScans,N_IR);
PPG_to_HR_IR = zeros(nScans,N_IR);
PPG_to_GS_IR = zeros(nScans,N_IR);
GS_to_Br_IR = zeros(nScans,N_IR);
GS_to_HR_IR = zeros(nScans,N_IR);
GS_to_PPG_IR = zeros(nScans,N_IR);

bounds_breathing = zeros(nScans,2,4);
bounds_HR = zeros(nScans,2,4);
bounds_PPG = zeros(nScans,2,4);
bounds_GS = zeros(nScans,2,4);

%%    Run loop across scans -------------------


tic
parfor c = 1 :  nScans
    
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    [Ts,Fs,TR, trig,PPGlocs, HR, resp, DVARS, cardiac, GS, FD, RF, BR, movRegr]  =  load_scan(subject,task,baseDir);
    
    HRmean(c)  = mean(HR);     HRstd(c) = std(HR);    FDmean(c) = mean(FD);
    GSmean(c) = mean(GS);     GSstd(c) = std(GS);
    
    volDel=0;
    time = 0:Ts:(length(trig)-1)*Ts;
    time_10 = 0:0.1:time(end);
    ind_BOLD=find(trig==1);
    trig(ind_BOLD(1:volDel))=0;
    ind_BOLD=find(trig==1);
    timeMR=time(trig==1); timeMR = timeMR +TR/2;
    NV = length(timeMR);
    GS = zscore(GS);
    
    
    %% Extract RV
    
    RV = zeros(size(resp));
    N = length(resp);    
    for i = 2:N-1
        ind_1 = i-6*Fs_10;   ind_1 = max(1,ind_1);
        ind_2 = i+6*Fs_10;   ind_2 = min(ind_2, N);
        RV(i) = std(resp(ind_1:ind_2));
    end
    RV(1) = RV(2); RV(end) = RV(end-1);    RV = RV(:);
    
    %% upsample GS
    
    t_scale = [0; timeMR(:); time(end)];
    x_scale = [GS(1); GS(:); GS(end)];
    GS_10 = interp1(t_scale, x_scale, time_10);
    
    time_10 = time_10(:);    GS_10 = GS_10(:);     resp = resp(:);
    %     RF = diff(resp); RF = [RF(1); RF]; RF = RF.^2;    
    
    %%  ----------------------------
    
    
    output = [resp, HR, cardiac, GS_10];
    
    % Breathing
    input = resp;
    [CC_all, bounds_all] = estimate_CC(input, output, N_IR);
    Br_to_Br_IR(c,:) = CC_all(:,1);
    Br_to_HR_IR(c,:) = CC_all(:,2);
    Br_to_PPG_IR(c,:) = CC_all(:,3);
    Br_to_GS_IR(c,:) = CC_all(:,4);
    bounds_breathing(c,:,:) = bounds_all;
    
    % Heart rate
    input = HR;
    [CC_all, bounds_all] = estimate_CC(input, output, N_IR);
    HR_to_Br_IR(c,:) =  CC_all(:,1);
    HR_to_HR_IR(c,:) = CC_all(:,2);
    HR_to_PPG_IR(c,:) = CC_all(:,3);
    HR_to_GS_IR(c,:) = CC_all(:,4);
    bounds_HR(c,:,:) = bounds_all;
    
    % PPG
    input = cardiac;
    [CC_all, bounds_all] = estimate_CC(input, output, N_IR);
    PPG_to_Br_IR(c,:) = CC_all(:,1);
    PPG_to_HR_IR(c,:) = CC_all(:,2);
    PPG_to_PPG_IR(c,:) = CC_all(:,3);
    PPG_to_GS_IR(c,:) = CC_all(:,4);
    bounds_PPG(c,:,:) = bounds_all;
    
    % GS
    input = GS_10;
    [CC_all, bounds_all] = estimate_CC(input, output, N_IR);
    GS_to_Br_IR(c,:) = CC_all(:,1);
    GS_to_HR_IR(c,:) = CC_all(:,2);
    GS_to_PPG_IR(c,:) = CC_all(:,3);
    GS_to_GS_IR(c,:) = CC_all(:,4);
    bounds_GS(c,:,:) = bounds_all;
    
    % RF
    input = RF;
    output = [RF, HR, cardiac, GS_10];
    [CC_all, bounds_all] = estimate_CC(input, output, N_IR);
    RF_to_RF_IR(c,:) = CC_all(:,1);
    RF_to_HR_IR(c,:) = CC_all(:,2);
    RF_to_PPG_IR(c,:) = CC_all(:,3);
    RF_to_GS_IR(c,:) = CC_all(:,4);
    bounds_RF(c,:,:) = bounds_all;

        
    % RV
    input = RV;
    output = [RV, HR, cardiac, GS_10];
    [CC_all, bounds_all] = estimate_CC(input, output, N_IR);
    RV_to_RV_IR(c,:) = CC_all(:,1);
    RV_to_HR_IR(c,:) = CC_all(:,2);
    RV_to_PPG_IR(c,:) = CC_all(:,3);
    RV_to_GS_IR(c,:) = CC_all(:,4);
    bounds_RV(c,:,:) = bounds_all;
    
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60), load chirp,  sound(y,Fs)

save([save_path,path_filename])


%%  Plot timeseris for the case with highest number of peaks in HR due to DVARS


nPeaks_HR_all_subj = reshape(nPeaks_HR_all,[4,nSubj])';
HRmean_subj = reshape(HRmean,[4,nSubj])';

fprintf('ICC for # of HR peaks: %3.2f;     \n',ICC(nPeaks_DVARS_before_HR_all_subj,'1-k'))


%% Average within subject

Br_to_Br_IR_sbj = IR_to_IRsbj(Br_to_Br_IR);
Br_to_HR_IR_sbj = IR_to_IRsbj(Br_to_HR_IR);
Br_to_PPG_IR_sbj = IR_to_IRsbj(Br_to_PPG_IR);
Br_to_GS_IR_sbj = IR_to_IRsbj(Br_to_GS_IR);

HR_to_Br_IR_sbj = IR_to_IRsbj(HR_to_Br_IR);
HR_to_HR_IR_sbj = IR_to_IRsbj(HR_to_HR_IR);
HR_to_PPG_IR_sbj = IR_to_IRsbj(HR_to_PPG_IR);
HR_to_GS_IR_sbj = IR_to_IRsbj(HR_to_GS_IR);

PPG_to_Br_IR_sbj = IR_to_IRsbj(PPG_to_Br_IR);
PPG_to_HR_IR_sbj = IR_to_IRsbj(PPG_to_HR_IR);
PPG_to_PPG_IR_sbj = IR_to_IRsbj(PPG_to_PPG_IR);
PPG_to_GS_IR_sbj = IR_to_IRsbj(PPG_to_GS_IR);

GS_to_Br_IR_sbj = IR_to_IRsbj(GS_to_Br_IR);
GS_to_HR_IR_sbj = IR_to_IRsbj(GS_to_HR_IR);
GS_to_PPG_IR_sbj = IR_to_IRsbj(GS_to_PPG_IR);
GS_to_GS_IR_sbj = IR_to_IRsbj(GS_to_GS_IR);

RF_to_RF_IR_sbj = IR_to_IRsbj(RF_to_RF_IR);
RF_to_HR_IR_sbj = IR_to_IRsbj(RF_to_HR_IR);
RF_to_PPG_IR_sbj = IR_to_IRsbj(RF_to_PPG_IR);
RF_to_GS_IR_sbj = IR_to_IRsbj(RF_to_GS_IR);

RV_to_RV_IR_sbj = IR_to_IRsbj(RV_to_RV_IR);
RV_to_HR_IR_sbj = IR_to_IRsbj(RV_to_HR_IR);
RV_to_PPG_IR_sbj = IR_to_IRsbj(RV_to_PPG_IR);
RV_to_GS_IR_sbj = IR_to_IRsbj(RV_to_GS_IR);


%%  ---------------------------------------------

nScans = 400;
t_IR = -50:0.1:50;
line_origin_linewidth = 2;

figure

% Breathing
ax(1) = subplot(4,4,1);
bounds = mean(bounds_breathing(:,1,1));
IR_all = Br_to_Br_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'g',0.9) , hold on, grid on
% title('Autocorrelation in breathing')

ax(5) = subplot(4,4,5);
bounds = mean(bounds_breathing(:,1,2));
IR_all = Br_to_HR_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Breathing --> Heart rate')

ax(9) = subplot(4,4,9);
bounds = mean(bounds_breathing(:,1,3));
IR_all = Br_to_PPG_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Breathing --> PPG-Amp')

ax(13) = subplot(4,4,13);
bounds = mean(bounds_breathing(:,1,4));
IR_all = Br_to_GS_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Breathing --> GS')

% Heart rate
ax(2) = subplot(4,4,2);
bounds = mean(bounds_HR(:,1,1));
IR_all = HR_to_Br_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Heart rate --> Breathing')

ax(6) = subplot(4,4,6);
bounds = mean(bounds_HR(:,1,2));
IR_all = HR_to_HR_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'g',0.9) , hold on, grid on
% title('Autocorrelation in heart rate')

ax(10) = subplot(4,4,10);
bounds = mean(bounds_HR(:,1,3));
IR_all = HR_to_PPG_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Heart rate --> PPG-Amp')

ax(14) = subplot(4,4,14);
bounds = mean(bounds_HR(:,1,4));
IR_all = HR_to_GS_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Heart rate --> GS')

% PPG
ax(3) = subplot(4,4,3);
bounds = mean(bounds_PPG(:,1,1));
IR_all = PPG_to_Br_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('PPG-Amp --> Breathing')

ax(7) = subplot(4,4,7);
bounds = mean(bounds_PPG(:,1,2));
IR_all = PPG_to_HR_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('PPG-Amp --> Heart rate')

ax(11) = subplot(4,4,11);
bounds = mean(bounds_PPG(:,1,3));
IR_all = PPG_to_PPG_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'g',0.9) , hold on, grid on
% title('Autocorrelation in PPG-Amp')

ax(15) = subplot(4,4,15);
bounds = mean(bounds_PPG(:,1,4));
IR_all = PPG_to_GS_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('PPG-Amp --> GS')

% Global signal
ax(4) = subplot(4,4,4);
bounds = mean(bounds_GS(:,1,1));
IR_all = GS_to_Br_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('GS --> Breathing')

ax(8) = subplot(4,4,8);
bounds = mean(bounds_GS(:,1,2));
IR_all = GS_to_HR_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('GS --> Heart rate')

ax(12) =  subplot(4,4,12);
bounds = mean(bounds_GS(:,1,3));
IR_all = GS_to_PPG_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('GS --> PPG-Amp')

ax(16) = subplot(4,4,16);
bounds = mean(bounds_GS(:,1,4));
IR_all = GS_to_GS_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'g',0.9) , hold on, grid on
% title('Autocorrelation in GS')

linkaxes(ax, 'x')
xlim([-30 30])

for i = [13:16]
    subplot(ax(i))
    xlabel('Lag time (s)')
end

for i = [1,5,9,13]
    subplot(ax(i)); ylabel ('Cross-correlation')
end

for i = 1:4
    for j=1:4
        if i~=j
            k = (i-1)*4+j;
            subplot(ax(k));
            ylim([-.3 .3])
        end
    end
end



%%  Replace Breathing with RV      ---------------------------------------------

nScans = 400;
t_IR = -50:0.1:50;
line_origin_linewidth = 2;

nSubj = 1;

figure

% Breathing
ax(1) = subplot(4,4,1);
bounds = mean(bounds_RV(:,1,1));
IR_all = RV_to_RV_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'g',0.9) , hold on, grid on
% title('Autocorrelation in breathing')

ax(5) = subplot(4,4,5);
bounds = mean(bounds_RV(:,1,2));
IR_all = RV_to_HR_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on


% title('Breathing --> Heart rate')

ax(9) = subplot(4,4,9);
bounds = mean(bounds_RV(:,1,3));
IR_all = RV_to_PPG_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Breathing --> PPG-Amp')

ax(13) = subplot(4,4,13);
bounds = mean(bounds_RV(:,1,4));
IR_all = RV_to_GS_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Breathing --> GS')


ax(6) = subplot(4,4,6);
bounds = mean(bounds_HR(:,1,2));
IR_all = HR_to_HR_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'g',0.9) , hold on, grid on
% title('Autocorrelation in heart rate')

ax(10) = subplot(4,4,10);
bounds = mean(bounds_HR(:,1,3));
IR_all = HR_to_PPG_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Heart rate --> PPG-Amp')

ax(14) = subplot(4,4,14);
bounds = mean(bounds_HR(:,1,4));
IR_all = HR_to_GS_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Heart rate --> GS')


ax(11) = subplot(4,4,11);
bounds = mean(bounds_PPG(:,1,3));
IR_all = PPG_to_PPG_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'g',0.9) , hold on, grid on
% title('Autocorrelation in PPG-Amp')

ax(15) = subplot(4,4,15);
bounds = mean(bounds_PPG(:,1,4));
IR_all = PPG_to_GS_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('PPG-Amp --> GS')


ax(16) = subplot(4,4,16);
bounds = mean(bounds_GS(:,1,4));
IR_all = GS_to_GS_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'g',0.9) , hold on, grid on
% title('Autocorrelation in GS')

linkaxes(ax, 'x')
xlim([-30 30])

for i = [13:16]
    subplot(ax(i))
    xlabel('Lag time (s)')
end

for i = [1,5,9,13]
    subplot(ax(i)); ylabel ('Cross-correlation')
end

% for i = [5,9,10,13,14,15]
%     subplot(ax(i));
%     ylim([-.45 .4])
% end




%%    RF and RV    ---------------------------------------------


nScans = 400;
t_IR = -50:0.1:50;
line_origin_linewidth = 2;

figure

% RF
ax(1) = subplot(4,2,1);
bounds = mean(bounds_RF(:,1,1));
IR_all = RF_to_RF_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'g',0.9) , hold on, grid on
% title('Autocorrelation in breathing')

ax(3) = subplot(4,2,3);
bounds = mean(bounds_RF(:,1,2));
IR_all = RF_to_HR_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on


% title('Breathing --> Heart rate')

ax(5) = subplot(4,2,5);
bounds = mean(bounds_RF(:,1,3));
IR_all = RF_to_PPG_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Breathing --> PPG-Amp')

ax(7) = subplot(4,2,7);
bounds = mean(bounds_RF(:,1,4));
IR_all = RF_to_GS_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Breathing --> GS')

% RV
ax(2) = subplot(4,2,2);
bounds = mean(bounds_RV(:,1,1));
IR_all = RV_to_RV_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'g',0.9) , hold on, grid on
% title('Heart rate --> Breathing')

ax(4) = subplot(4,2,4);
bounds = mean(bounds_RV(:,1,2));
IR_all = RV_to_HR_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Autocorrelation in heart rate')

ax(6) = subplot(4,2,6);
bounds = mean(bounds_RV(:,1,3));
IR_all = RV_to_PPG_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Heart rate --> PPG-Amp')

ax(8) = subplot(4,2,8);
bounds = mean(bounds_RV(:,1,4));
IR_all = RV_to_GS_IR_sbj;
IR = mean(IR_all); IRstd = std(IR_all/sqrt(nSubj));
xline(0,'LineWidth',line_origin_linewidth); hold on
yline(bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
yline(-bounds,'LineWidth',line_origin_linewidth,'linestyle',':'); hold on
shadedErrorBar(t_IR,IR,IRstd, 'r',0.9) , hold on, grid on
% title('Heart rate --> GS')


linkaxes(ax, 'x')
xlim([-20 40])

for i = [7:8]
    subplot(ax(i))
    xlabel('Lag time (s)')
end

for i = [1,3,5,7]
    subplot(ax(i)); ylabel ('Cross-correlation')
end

for i = 2:4
    for j=1:2
       
            k = (i-1)*2+j;
            subplot(ax(k));
            ylim([-.4 .4])
        end
    
end


%%  ---------------------------------------------


