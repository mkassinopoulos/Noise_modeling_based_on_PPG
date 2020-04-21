

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
HRmean_all = zeros(nScans,1);
HRstd_all = zeros(nScans,1);

r_CPM_all = zeros(nScans,nShifts,M_order);
r_CPM_Amp_all = zeros(nScans,nShifts,M_order);
r_RETR_all = zeros(nScans,nShifts,M_order);

tic
parfor c = 1 : nScans
    
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    filepath_MRacq=[baseDir,'/Physio/',subject,'_',task,'/phys.mat'];
    [Ts,Fs, PPGlocs, HR,  cardiac] =  load_scan_only_cardiac(subject,task,baseDir,filepath_MRacq);
    
    HRmean = mean(HR);
    
    voxel = cardiac;
    voxel(1:nDel) = []; voxel(end-nDel:end) = [];
    NV = length(voxel);
    voxel = filtfilt(filt_b,filt_a,voxel);
    N = length(cardiac); time = 0:Ts:(N-1)*Ts;
    
    block_length=round((NV)/kFold); ind_blocks=zeros(kFold,2);
    for i=1:kFold
        ind_blocks(i,1)=1+(i-1)*block_length;
        ind_blocks(i,2)=min(i*block_length,NV);  ind_blocks(kFold,2) = NV;
    end
    
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
    
    %%  ----------------------------------------------
    
    r_CPM_shift = zeros(nShifts,M_order);
    r_CPM_Amp_shift = zeros(nShifts,M_order);
    r_RETR_shift = zeros(nShifts,M_order);
    for  c_temp = 1:nShifts
        shift = temp_shift_scale(c_temp);
        
        CPM_regr = zeros(NV,M_order);
        CPM_Amp_regr = zeros(NV,M_order);
        RETR_regr = zeros(NV,M_order);
        for m = 1:M_order*2
            ind = 1:NV; ind = round(ind + nDel + shift*Fs);
            CPM_regr(:,m) = CPM_regr_all(ind,m);
            CPM_Amp_regr(:,m) = CPM_Amp_regr_all(ind,m);
            RETR_regr(:,m) = RETR_regr_all(ind,m);
        end
        
        regr_CPM = [CPM_regr, ones(NV,1)];
        regr_CPM_Amp = [CPM_Amp_regr, ones(NV,1)];
        RETR_regr = [RETR_regr, ones(NV,1)];
        
        %% CPM  ---------------------------------------------
        
        regr = regr_CPM;
        r_valid_k=zeros(M_order,kFold);
        for k=1:kFold
            ind_valid = ind_blocks(k,1):ind_blocks(k,2);
            ind_train = 1:NV; ind_train(ind_valid)=[];
            if kFold==1, ind_train=1:NV; end
            
            regr_train = regr(ind_train,:);  % regr_train(:,1:end-1) = detrend(regr_train(:,1:end-1),'linear');
            regr_valid = regr(ind_valid,:); % regr_valid(:,1:end-1) = detrend(regr_valid(:,1:end-1),'linear');
            voxel_train =voxel(ind_train);    voxel_valid = voxel(ind_valid);
            %             voxel_train = detrend(voxel(ind_train),'linear');    voxel_valid = detrend(voxel(ind_valid),'linear');
            
            for m = 1:M_order
                ind_regr = [1:2*m,size(regr,2)];
                B = regr_train(:,ind_regr)\voxel_train;     yPred = regr_valid(:,ind_regr)*B;
                r_valid_k(m,k) = corr(voxel_valid,yPred) ;
                
                if 0
                    timeVoxel = T_del:Ts:T_del+(NV-1)*Ts;
                    figure('Position', [ 72         868        2282         420])
                    plot(timeVoxel,voxel), hold on,       plot(timeVoxel,yPred)
                    %                     xlim([60 180])
                    xlabel('Time (s)'),             ylabel('Amplitude (a.u.)')
                    legend('Experimental data','Fitted data with DCPM')
                end
            end
        end
        r_CPM_shift(c_temp,:) = mean(r_valid_k,2);
        
        %% CPM_Amp ---------------------------------------------
        
        regr = regr_CPM_Amp;
        r_valid_k=zeros(M_order,kFold);
        for k=1:kFold
            ind_valid = ind_blocks(k,1):ind_blocks(k,2);
            ind_train = 1:NV; ind_train(ind_valid)=[];
            if kFold==1, ind_train=1:NV; end
            
            regr_train = regr(ind_train,:); % regr_train(:,1:end-1) = detrend(regr_train(:,1:end-1),'linear');
            regr_valid = regr(ind_valid,:); % regr_valid(:,1:end-1) = detrend(regr_valid(:,1:end-1),'linear');
            voxel_train =voxel(ind_train);    voxel_valid = voxel(ind_valid);
            %             voxel_train = detrend(voxel(ind_train),'linear');    voxel_valid = detrend(voxel(ind_valid),'linear');
            
            for m = 1:M_order
                ind_regr = [1:2*m,size(regr,2)];
                B = regr_train(:,ind_regr)\voxel_train;     yPred = regr_valid(:,ind_regr)*B;
                r_valid_k(m,k) = corr(voxel_valid,yPred);
                
                if 0
                    timeVoxel = T_del:Ts:T_del+(NV-1)*Ts;
                    figure('Position', [ 72         868        2282         420])
                    plot(timeVoxel,voxel), hold on,       plot(timeVoxel,yPred)
                    %                     xlim([60 180])
                    xlabel('Time (s)'),             ylabel('Amplitude (a.u.)')
                    legend('Experimental data','Fitted data with DCPM')
                end
                
            end
            
        end
        r_CPM_Amp_shift(c_temp,:) = mean(r_valid_k,2);
        
        %% RETROICOR  ---------------------------------------------
        
        regr = RETR_regr;
        r_valid_k=zeros(M_order,kFold);
        for k=1:kFold
            ind_valid = ind_blocks(k,1):ind_blocks(k,2);
            ind_train = 1:NV; ind_train(ind_valid)=[];
            if kFold==1, ind_train=1:NV; ind_valid = ind_train; end
            
            regr_train = regr(ind_train,:);  %  regr_train(:,1:end-1) = detrend(regr_train(:,1:end-1),'linear');
            regr_valid = regr(ind_valid,:); % regr_valid(:,1:end-1) = detrend(regr_valid(:,1:end-1),'linear');
            voxel_train =voxel(ind_train);    voxel_valid = voxel(ind_valid);
            %             voxel_train = detrend(voxel(ind_train),'linear');    voxel_valid = detrend(voxel(ind_valid),'linear');
            
            for m = 1:M_order
                ind_regr = [1:2*m,size(regr,2)];
                B = regr_train(:,ind_regr)\voxel_train;     yPred = regr_valid(:,ind_regr)*B;
                r_valid_k(m,k) = corr(voxel_valid,yPred);
                
                if 0
                    timeVoxel = T_del:Ts:T_del+(NV-1)*Ts;
                    figure('Position', [ 72         868        2282         420])
                    plot(timeVoxel,voxel), hold on,       plot(timeVoxel,yPred)
                    %                      xlim([60 180])
                    xlabel('Time (s)'),             ylabel('Amplitude (a.u.)')
                    legend('Experimental data','Fitted data with RETROICOR')
                end
            end
            
        end
        r_RETR_shift(c_temp,:) = mean(r_valid_k,2);
        
        
    end
    r_CPM_all(c,:,:) = r_CPM_shift;
    r_CPM_Amp_all(c,:,:) = r_CPM_Amp_shift;
    r_RETR_all(c,:,:) = r_RETR_shift;
end

save([save_path,'workspace'], 'r_CPM_all','r_CPM_Amp_all','r_RETR_all','subject_list')

%%  ---------------------

figure

x = squeeze(r_RETR_all(:,:,8));
NP = length(temp_shift_scale);
x_sbj = zeros(nSubj,NP);
for i = 1:nSubj
    ind = (i-1)*4 + (1:4);
    tmp = x(ind,:);
    x_sbj(i,:) = mean(tmp,1);
end
RETR_sbj = x_sbj;
r_mean_RETR = mean(x_sbj);
r_std_RETR = std(x_sbj)/sqrt(nSubj);



x = squeeze(r_CPM_all(:,:,8));
NP = length(temp_shift_scale);
x_sbj = zeros(nSubj,NP);
for i = 1:nSubj
    ind = (i-1)*4 + (1:4);
    tmp = x(ind,:);
    x_sbj(i,:) = mean(tmp,1);
end
CPM_CA_sbj = x_sbj;
r_mean_CPM_CA = mean(x_sbj);
r_std_CPM_CA = std(x_sbj)/sqrt(nSubj);



x = squeeze(r_CPM_Amp_all(:,:,8));
NP = length(temp_shift_scale);
x_sbj = zeros(nSubj,NP);
for i = 1:nSubj
    ind = (i-1)*4 + (1:4);
    tmp = x(ind,:);
    x_sbj(i,:) = mean(tmp,1);
end
CPM_VA_sbj = x_sbj;
r_mean_CPM_VA = mean(x_sbj);
r_std_CPM_VA = std(x_sbj)/sqrt(nSubj);


figure
shadedErrorBar(-temp_shift_scale,r_mean_RETR,r_std_RETR, 'r') , hold on
shadedErrorBar(-temp_shift_scale,r_mean_CPM_CA,r_std_CPM_CA, 'k',0.9) , hold on
shadedErrorBar(-temp_shift_scale,r_mean_CPM_VA,r_std_CPM_VA, 'g',0.9) , hold on
grid on
xlim([-2 2])
xlabel('Lag time (s)'), ylabel('Mean correlation')

%%  -----------------

x1 = RETR_sbj(:,32);
x2 = CPM_CA_sbj(:,33);
x3 = CPM_VA_sbj(:,36);



mean(x1)
mean(x2)
mean(x3)

[ttest_h ttest_p] = ttest(x1,x2)

%%  ------------------------------

dx_sbj = x2-x1;
dx_relat = (x2-x1)*100./x1;


HRstd_sbj = reshape(HRstd_all,[4 100])'; HRstd_sbj = mean(HRstd_sbj,2);


figure

x = HRstd_sbj; 

y = dx_sbj;
y = dx_relat;
scatter(x,y,50,'k','LineWidth',2)
r = corr(x,y)

model = fitlm(x,y)
pVal = model.Coefficients.pValue(2);
ax=gca; ax.FontWeight='bold'; hold on
Par = model.Coefficients.Estimate; P12 = Par(1)+Par(2)*1;
x_scal = (min(x)-0.2:0.0001:max(x)+0.2);
y_scal = Par(1) + Par(2)*x_scal;
plot(x_scal,y_scal,'r','LineWidth',3)

xlabel('Heart rate variability (bpm)')
ylabel('?r(RETR-->CPM_{CA})')
ylabel('% improvement')

grid on

xlim([1 10.5])
% ylim([-0.013 0.022])


%%  ------------

if 0
    figure('position',[800 400 1500 900])
    
    data = [x3,x2, x1];
    
    
    subplot(2,2,pl)
    boxplot(data), [xl] = get(gca,'XLim');   [yl] = get(gca,'YLim');
    
    bx_labels = {'RETROICOR','DCPM_{CA}','DCPM_{VA}'};
    boxplot(data,'Labels',bx_labels)
    
    ax = gca;    ax.YGrid = 'on';    ax.YGrid;    ax.GridAlpha=0.3;
    ax.GridLineStyle='--';    ax.FontSize=14;
    title(titleTxt),    xlabel('Model'),     ylabel('Correlation (%)')
    ax.XTickLabelRotation = 50
    ax.FontWeight = 'bold'
end






