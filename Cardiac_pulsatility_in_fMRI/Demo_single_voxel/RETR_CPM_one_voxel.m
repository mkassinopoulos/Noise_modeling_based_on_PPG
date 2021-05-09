%% Load sample data.
% Physio data for HCP subject S515014 (Rest1LR) and
%  timeseries of a voxel 

clear, clc, close all

load('sample_phys_voxel_S105014_R1LR.mat')   % HCP subject S515014 (Rest1LR)

addpath('..')

volDel=20; volDel_from_end=10;

TR = 0.72;  Ts = 1/Fs;
M_order = 6; % Order of fourier series for models

HRmean = mean(HR_10);
N = length(cardiac);

ind_BOLD=find(trig==1);
trig(ind_BOLD(1:volDel))=0;    trig(ind_BOLD(end-volDel_from_end+1:end)) = 0;
ind_BOLD=find(trig==1);

time = 0:Ts:(length(trig)-1)*Ts;
timeMR=time(trig==1);
NV = length(timeMR);

voxel = zscore(voxel(1:end-volDel_from_end));

nPeaks = length(PPGlocs); 
ind_PPGlocs = zeros(nPeaks,1);
for i = 1:nPeaks
    t = PPGlocs(i);
    [val loc] = min(abs(time-t));
    ind_PPGlocs(i) = loc;
end

c_green = [0.47,0.67,0.19];
c_orange = [0.85,0.33,0.10];
c_gray = [0.80,0.80,0.80];
c_black = [0.00,0.00,0.00];

time_10 = 0:0.1:(length(HR_10)-1)*0.1;
plot(time_10, HR_10)
xlabel('Time (s)'), ylabel('Heart rate (bpm)')

%% S1a: RETROICOR

[regr_card, Phi] = RETR_Card_regressors_v2(time,PPGlocs,M_order);
shift = 0;
regr_card_MR = regr_card(ind_BOLD+round(shift*Fs),:);
regr_card_high_Fs = [regr_card,ones(length(time),1)];
regr_all_MR = [regr_card_MR, ones(NV,1)];


B = regr_all_MR\voxel ;   
yPred =regr_all_MR*B;
yPred_high_Fs = regr_card_high_Fs*B;
r  = corr(voxel,yPred)

figure('position',[ 265         547        1431         420])
plot(voxel), hold on
plot(yPred)
xlabel('Time (s)'), ylabel('Voxel timeseries (a.u.)')
title(sprintf('RETROICOR (Correlation r: %3.2f) ',r))

figure('position',[  435         558        1167         420])
plot(time,regr_card)
xlim([400 404])
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')
title('Nuisance regressors at high sampling rate (RETROICOR)')

Phi_scale = 0:0.01:2*pi;
for i = 1: M_order
    RETR_basis(:,(i-1)*2+1) = cos(i*Phi_scale);
    RETR_basis(:,(i-1)*2+2) = sin(i*Phi_scale);
end
CPW_RETR = RETR_basis*B(1:2*M_order);
CPW_RETR = CPW_RETR/max(abs(CPW_RETR));

figure('Position',[ 1114         554         467         322])
plot(Phi_scale,CPW_RETR,'Color',c_orange,'LineWidth',3)
grid on
xlabel('Cardiac phase (rad)')
ylabel('Amplitude (a.u.)')
title('Cardiac pulsatility waveform (CPW)')
xlim([0 max(Phi_scale)])
ylim([-1.2 1.2])
set(gca,'XAxisLocation','origin'), box off
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})


%%  S1b: Plot regr - Card-RETROICOR
% Demonstration of Card-RETROICOR as shown in Suppl. Fig. 1 (Kassinopoulos & Mitsis, 2021)

x1 = 400; x2 = 404;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontWeight','normal')

close all

fig = figure('Position',[1335          63         585         933]);
left_color = c_green;
right_color =  c_orange;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

ax1 = subplot(7,1,1);
yyaxis left
plot(time,cardiac,'Color', c_green,'LineWidth',2), hold on
h1 = plot(PPGlocs, cardiac(ind_PPGlocs),'o','Color',c_green);
set(h1, 'markerfacecolor', get(h1, 'color')); 
xticks([400:404]), ylim([-3.5 3.5 ])

yyaxis right
plot(time, Phi,'.')
ylim([-2.5 2*pi+2.5 ])
yticks([0 pi 2*pi])
yticklabels({'0','\pi','2\pi'})

ax2 = subplot(7,1,2);
plot(time, regr_card(:,1),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,1),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
% set(h1, 'markerfacecolor', get(h1, 'color')); 
xticks([400:404])

ax3 = subplot(7,1,3);
plot(time, regr_card(:,2),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,2),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
xticks([400:404])

ax4 = subplot(7,1,4);
plot(time, regr_card(:,3),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,3),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
xticks([400:404])

ax5 = subplot(7,1,5);
plot(time, regr_card(:,4),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,4),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
xticks([400:404])

ax6 = subplot(7,1,6);
plot(time, yPred_high_Fs,'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, voxel,'-*','Color','k','LineWidth',2,'MarkerSize',14); hold on
plot(timeMR, yPred,'-*','Color',right_color,'LineWidth',2,'MarkerSize',14); hold on

xticks([400:404])

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6], 'x')
xlim([x1 x2])

ax7 = subplot(7,1,7);
plot(timeMR, voxel,'-','Color','k','LineWidth',2); hold on
plot(timeMR, yPred,'-','Color',right_color,'LineWidth',2); hold on
xlim([350 450])



%% S2a: CPM_ca

u = zeros(size(time)); uA = u;
nPeaks = length(PPGlocs);
for i = 1:nPeaks
    t = PPGlocs(i);
    [val loc] = min(abs(time-t));
    u(loc) = 1;      uA(loc) = cardiac(loc);
end

memory = 60/HRmean;
CPM_IR = func_CPM_cos(Ts, memory, M_order);

CPM_regr_all = zeros(N,M_order);
CPM_Amp_regr_all = zeros(N,M_order);
for m = 1:M_order*2
    x = conv(u,CPM_IR(:,m));    CPM_regr_all(:,m) = x(1:N);
    x = conv(uA,CPM_IR(:,m));  CPM_Amp_regr_all(:,m) = x(1:N);
end

shift = 0;
regr_card = CPM_regr_all;
regr_card_MR = regr_card(ind_BOLD+round(shift*Fs),:);
regr_all_MR = [regr_card_MR, ones(NV,1)];
regr_card_high_Fs = [regr_card,ones(length(time),1)];

B = regr_all_MR\voxel ;     
yPred =regr_all_MR*B;
yPred_high_Fs = regr_card_high_Fs*B;
r  = corr(voxel,yPred)

figure('position',[ 265         547        1431         420])
plot(voxel), hold on
plot(yPred)
title(sprintf('CPM_{CA} (Correlation r: %3.2f) ',r))

figure('Position',[ 1114         554         467         322])
CPW = CPM_IR*B(1:2*M_order); CPW = CPW/max(abs(CPW));
t_IR = 0:Ts:(length(CPW)-1)*Ts;
plot(t_IR,CPW,'Color',c_orange,'LineWidth',3)
grid on
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
title('Cardiac pulsatility waveform (CPW)')
xlim([0 max(t_IR)])
ylim([-1.2 1.2])
set(gca,'XAxisLocation','origin'), box off

%%  S2b: Plot regr for CPM_CA
% Demonstration of CPM_CA as shown in Suppl. Fig. 1 (Kassinopoulos & Mitsis, 2021)

x1 = 400; x2 = x1+4;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontWeight','normal')

close all

fig = figure('Position',[1335          63         585         933]);
left_color = c_green;
right_color =  c_orange;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

ax1 = subplot(7,1,1);
yyaxis left
plot(time,cardiac,'Color', c_green,'LineWidth',2), hold on
h1 = plot(PPGlocs, cardiac(ind_PPGlocs),'o','Color',c_green);
set(h1, 'markerfacecolor', get(h1, 'color')); 
xticks([x1:x2]), ylim([-3.5 3.5 ])

yyaxis right
plot(time, u,'-')
xticks([x1:x2]), ylim([-3.5 3.5 ])

ax2 = subplot(7,1,2);
plot(time, regr_card(:,1),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,1),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
% set(h1, 'markerfacecolor', get(h1, 'color')); 
xticks([x1:x2])

ax3 = subplot(7,1,3);
plot(time, regr_card(:,2),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,2),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
xticks([x1:x2])

ax4 = subplot(7,1,4);
plot(time, regr_card(:,3),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,3),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
xticks([x1:x2])

ax5 = subplot(7,1,5);
plot(time, regr_card(:,4),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,4),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
xticks([x1:x2])

ax6 = subplot(7,1,6);
plot(time, yPred_high_Fs,'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, voxel,'-*','Color','k','LineWidth',2,'MarkerSize',14); hold on
plot(timeMR, yPred,'-*','Color',right_color,'LineWidth',2,'MarkerSize',14); hold on
xticks([x1:x2])

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6], 'x')
xlim([x1 x2])

ax7 = subplot(7,1,7);
plot(timeMR, voxel,'-','Color','k','LineWidth',2); hold on
plot(timeMR, yPred,'-','Color',right_color,'LineWidth',2); hold on
xlim([350 450])


%% S3a: CPM_va

shift = 0;
regr_card = CPM_Amp_regr_all;
regr_card_MR = regr_card(ind_BOLD+round(shift*Fs),:);
regr_all_MR = [regr_card_MR, ones(NV,1)];
regr_card_high_Fs = [regr_card,ones(length(time),1)];

B = regr_all_MR\voxel ;     
yPred =regr_all_MR*B;
yPred_high_Fs = regr_card_high_Fs*B;
r  = corr(voxel,yPred)

figure('position',[ 265         547        1431         420])
plot(voxel), hold on
plot(yPred)
title(sprintf('CPM_{VA} (Correlation r: %3.2f) ',r))

figure('Position',[ 1114         554         467         322])
CPW = CPM_IR*B(1:2*M_order); CPW = CPW/max(abs(CPW));
t_IR = 0:Ts:(length(CPW)-1)*Ts;
plot(t_IR,CPW,'Color',c_orange,'LineWidth',3)
grid on
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
title('Cardiac pulsatility waveform (CPW)')
xlim([0 max(t_IR)])
set(gca,'XAxisLocation','origin'), box off
ylim([-1.2 1.2])

%%  S3b: Plot regr for CPM_VA
% Demonstration of CPM_VA-CPM_VA as shown in Suppl. Fig. 1 (Kassinopoulos & Mitsis, 2021)

c_green = [0.47,0.67,0.19];
c_orange = [0.85,0.33,0.10];
c_gray = [0.80,0.80,0.80];
c_black = [0.00,0.00,0.00];

x1 = 400; x2 = 404;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',14)
set(0,'defaultAxesFontWeight','normal')

close all

fig = figure('Position',[1335          63         585         933]);
left_color = c_green;
right_color =  c_orange;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

ax1 = subplot(7,1,1);
yyaxis left
plot(time,cardiac,'Color', c_green,'LineWidth',2), hold on
h1 = plot(PPGlocs, cardiac(ind_PPGlocs),'o','Color',c_green);
set(h1, 'markerfacecolor', get(h1, 'color')); 
xticks([400:404]), ylim([-3.5 3.5 ])

yyaxis right
plot(time, uA,'-')
ylim([-3.5 3.5 ])

ax2 = subplot(7,1,2);
plot(time, regr_card(:,1),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,1),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
% set(h1, 'markerfacecolor', get(h1, 'color')); 
xticks([400:404])

ax3 = subplot(7,1,3);
plot(time, regr_card(:,2),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,2),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
xticks([400:404])

ax4 = subplot(7,1,4);
plot(time, regr_card(:,3),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,3),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
xticks([400:404])

ax5 = subplot(7,1,5);
plot(time, regr_card(:,4),'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, regr_card_MR(:,4),'-*','Color',right_color,'LineWidth',1,'MarkerSize',14);
xticks([400:404])

ax6 = subplot(7,1,6);
plot(time, yPred_high_Fs,'Color',c_gray,'LineWidth',2); hold on
plot(timeMR, voxel,'-*','Color','k','LineWidth',2,'MarkerSize',14); hold on
plot(timeMR, yPred,'-*','Color',right_color,'LineWidth',2,'MarkerSize',14); hold on
xticks([400:404])

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6], 'x')
xlim([x1 x2])

ax7 = subplot(7,1,7);
plot(timeMR, voxel,'-','Color','k','LineWidth',2); hold on
plot(timeMR, yPred,'-','Color',right_color,'LineWidth',2); hold on
xlim([350 450])











