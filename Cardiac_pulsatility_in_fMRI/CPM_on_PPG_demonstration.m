

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

save_path = ['Export/2020_02_01/CPM_on_PPG/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

Fs = 400; Ts = 1/Fs;  HPF_f = 0.008;  [filt_b,filt_a] = butter(2,HPF_f*2*Ts,'high');

kFold=3;
temp_shift_scale = -3:0.1:3; nShifts = length(temp_shift_scale);

%%    Run loop across scans -------------------

T_del = 5;   % delete first 5 seconds from output
nDel = T_del*Fs;
M_order = 8;

c = 74*4 +1;

s = ceil(c/4);        run = c - (s-1)*4;
subject = char(subject_list(s,:));         task = char(task_list(run));
fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)

filepath_MRacq=[baseDir,'/Physio/',subject,'_',task,'/phys.mat'];
[Ts,Fs, PPGlocs, HR,  cardiac] =  load_scan_only_cardiac(subject,task,baseDir,filepath_MRacq);

HRmean = mean(HR);
fprintf('HR: %3.1f+-%3.1f    \n',mean(HR),std(HR))

voxel = cardiac;
voxel(1:nDel) = []; voxel(end-nDel:end) = [];
NV = length(voxel);
voxel = filtfilt(filt_b,filt_a,voxel);
N = length(cardiac); time = 0:Ts:(N-1)*Ts;
time_10 = 0: 0.1 : time(end);  time_10 = time_10(:);

%%  ---------------------------------------------

memory = 60/HRmean;
CPM_IR = func_CPM_cos(Ts, memory, M_order);

u = zeros(size(time)); uA = u;
nPeaks = length(PPGlocs);
for i = 1:nPeaks
    t = PPGlocs(i);
    [val loc] = min(abs(time-t));
    u(loc) = 1;
    uA(loc) = cardiac(loc);
end

CPM_regr_all = zeros(N,M_order);
CPM_Amp_regr_all = zeros(N,M_order);
for m = 1:M_order*2
    u_conv = conv(u,CPM_IR(:,m)); u_conv = u_conv(1:N);
    CPM_regr_all(:,m) = u_conv(:);
    x = conv(uA,CPM_IR(:,m)); x = x(1:N);
    CPM_Amp_regr_all(:,m) = x(:);
end
CPM_regr_all = filtfilt(filt_b,filt_a,CPM_regr_all);
CPM_Amp_regr_all = filtfilt(filt_b,filt_a,CPM_Amp_regr_all);

RETR_regr_all = RETR_Card_regressors_v2(time,PPGlocs,M_order);
RETR_regr_all = filtfilt(filt_b,filt_a,RETR_regr_all);

%  ----------------------------------------------

shift_CPM = 0.5;
shift_CPM_VA = 0.5;
shift_RETR = 0;

ind = 1:NV; ind = round(ind + nDel + shift_CPM*Fs);
CPM_regr = CPM_regr_all(ind,:);

ind = 1:NV; ind = round(ind + nDel + shift_CPM_VA*Fs);
CPM_Amp_regr= CPM_Amp_regr_all(ind,:);

ind = 1:NV; ind = round(ind + nDel + shift_RETR*Fs);
RETR_regr= RETR_regr_all(ind,:);

regr_CPM = [CPM_regr, ones(NV,1)];
regr_CPM_Amp = [CPM_Amp_regr, ones(NV,1)];
RETR_regr = [RETR_regr, ones(NV,1)];

%% RETR  ---------------------------------------------

color_o = ['1.00,0.41,0.16'];

figure

timeVoxel = T_del:Ts:T_del+(NV-1)*Ts;

ax1 = subplot(4,1,1);
plot(time_10, HR,'k')
ylabel('Heart rate (bpm)')
grid

ax2 = subplot(4,1,2);
B = RETR_regr\voxel; yPred_RETR = RETR_regr*B;
plot(timeVoxel, voxel,'k'), hold on
plot(timeVoxel, yPred_RETR,'color', color_o)
ylabel('PPG (a.u.)')
grid
corr(voxel, yPred_RETR)

ax3 = subplot(4,1,3);
B = regr_CPM\voxel; yPred_CPM = regr_CPM*B;
plot(timeVoxel, voxel,'k'), hold on
plot(timeVoxel, yPred_CPM,'color', color_o)
ylabel('PPG (a.u.)')
grid
corr(voxel, yPred_CPM)


ax4 = subplot(4,1,4);
B = regr_CPM_Amp\voxel; yPred_CPM_VA = regr_CPM_Amp*B;
plot(timeVoxel, voxel,'k'), hold on
plot(timeVoxel, yPred_CPM_VA,'color', color_o)
ylabel('PPG (a.u.)')
%     title('CPM_{VA}')
grid
corr(voxel, yPred_CPM_VA)

linkaxes([ax1, ax2, ax3, ax4],'x')
xlabel('Time (s)')

xlim([ 400 430])

%%  ------------






















