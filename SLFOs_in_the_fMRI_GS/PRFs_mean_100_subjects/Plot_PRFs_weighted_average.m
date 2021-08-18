%% Plot PRFs with weighted average (no standard error) --------------------

clc, clear

flag_shift_ePPG = 1;   % set to 1 to use PARFs obtained with PPG-Amp shifted backward 5 seconds

if flag_shift_ePPG == 1
    load('PRFs_ePPG_shift_m5/PRFs_100_subjects.mat')
else
    load('PRFs_ePPG_shift_0/PRFs_100_subjects.mat')
end

nScans = 400;
N_win = size(CRF_all,1);
t_win = 0:0.1:(N_win-1)*0.1;

CRF_mean = zeros(N_win,1);
PARF_mean = zeros(N_win,1);
RRF_mean = zeros(N_win,1);
for c = 1:nScans
    x = CRF_all(:,c)*r_all(c,2);
    CRF_mean = CRF_mean + x/nScans;
    
    x = PARF_all(:,c)*r_all(c,4);
    PARF_mean = PARF_mean + x/nScans;
    
    x = RRF_all(:,c)*r_all(c,3);
    RRF_mean = RRF_mean + x/nScans;
end


figure('position', [  763   375   555   417])
ax = plot(t_win,[CRF_mean, PARF_mean, RRF_mean],'linewidth',3); grid on
ax = gca;
ax.GridLineStyle = '--';
ax.Box = 'off';
ax.XAxisLocation = 'origin';
legend('CRF','PARF','RRF');  legend boxoff
xlim([0 50]), ylim([-0.5 0.5])
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')


%% Plot PRFs with weighted average (and standard error) --------------------


CRF_all_subj = zeros(N_win,nScans/4);
PARF_all_subj = zeros(N_win,nScans/4);
RRF_all_subj = zeros(N_win,nScans/4);
r_card_subj = zeros(nScans/4,1);
r_brf_subj = zeros(nScans/4,1);
r_rrf_subj = zeros(nScans/4,1);
for s = 1:100
    ind = (1:4) + (s-1)*4;
    tmp = CRF_all(:,ind);     CRF_all_subj(:,s) = mean(tmp,2);
    tmp = PARF_all(:,ind);     PARF_all_subj(:,s) = mean(tmp,2);
    tmp = RRF_all(:,ind);     RRF_all_subj(:,s) = mean(tmp,2);
    
    r_card_subj(s) = mean(r_all(ind,2));
    r_brf_subj(s) = mean(r_all(ind,4));
    r_rrf_subj(s) = mean(r_all(ind,3));
end
CRF_std = std(CRF_all_subj', r_card_subj); CRF_std = CRF_std/sqrt(100);
BRF_std = std(PARF_all_subj',r_brf_subj); BRF_std = BRF_std/sqrt(100);
RRF_std = std(RRF_all_subj',r_rrf_subj); RRF_std = RRF_std/sqrt(100);


figure('position', [  763   375   555   417])

shadedErrorBar(t_win,CRF_mean,CRF_std, 'r',0.9) , hold on, grid on
shadedErrorBar(t_win,PARF_mean,BRF_std, 'k',0.9) , hold on, grid on
shadedErrorBar(t_win,RRF_mean,RRF_std, 'g',0.9) , hold on, grid on

grid on
ax = gca;
ax.GridLineStyle = '--';
ax.Box = 'off';
ax.XAxisLocation = 'origin';
xlim([0 50]), ylim([-0.5 0.5])
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')



