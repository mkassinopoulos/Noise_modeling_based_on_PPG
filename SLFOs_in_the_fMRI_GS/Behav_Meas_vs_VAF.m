


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

load('E:\CloudStation\HCP\Behavioral_Measurements\Beh_Meas_num.mat')

indSubj = zeros(nSubj,1);
for s = 1 :nSubj
    ind = find( BehMeas_num{:,1} == str2double(subject_list(s)));
    indSubj(s) = ind;
end    
BehMeas_num = BehMeas_num(indSubj,:);

%%    Mean of hematocrit

x = BehMeas_num{:,19:20};

hematocrit = zeros(nSubj,1);
for s = 1:nSubj
    tmp = x(s,:);    tmpMean = mean(tmp);
    if tmpMean~=NaN
        hematocrit(s) = tmpMean;
    else
        hematocrit(s) = max(tmp);
    end
end 


%%    Mean of blood pressure

x = BehMeas_num{:,21:22};

PressurePulse = zeros(nSubj,1);
for s = 1:nSubj
    tmp = x(s,:);    tmpMean = mean(tmp);
    if tmpMean~=NaN
        PressurePulse(s) = tmp(:,1)-tmp(:,2);
    else
        PressurePulse(s) = NaN;
    end
end 



%%   Confounds (HR, RV, motion)
Fs_10 = 10; Ts_10 = 0.1;


HRmean = zeros(nScans,1);
HRstd = zeros(nScans,1);
FDmean = zeros(nScans,1);
GSmean = zeros(nScans,1);
GSstd = zeros(nScans,1);
RVmean = zeros(nScans,1);
RVstd = zeros(nScans,1);
PPGmean = zeros(nScans,1);
PPGstd = zeros(nScans,1);

parfor c = 1   : nScans
       s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    [Ts,Fs,TR, trig,PPGlocs, HR, resp_10, DVARS, cardiac, GS, FD, ~, ~, movRegr]  =  load_scan(subject,task,baseDir);
    
    GS(1:20) = [];    
    HRmean(c)  = mean(HR);     HRstd(c) = std(HR);    FDmean(c) = mean(FD);
    GSmean(c) = mean(GS);     GSstd(c) = std(GS);
    PPGmean(c) = mean(cardiac);     PPGstd(c) = std(cardiac);
            
    %% Extract RV
    
    resp_10 = zscore(resp_10);
    RV = zeros(size(resp_10));
    N = length(resp_10);
    for i = 2:N-1
        ind_1 = i-6*Fs_10;   ind_1 = max(1,ind_1);
        ind_2 = i+6*Fs_10;   ind_2 = min(ind_2, N);
        RV(i) = std(resp_10(ind_1:ind_2));
    end
    RV(1) = RV(2); RV(end) = RV(end-1);    RV = RV(:);    
    RVmean(c) = mean(RV); RVstd(c) = std(RV);
    
end

save('Stats_on_Phys_100_subjects','HRmean','HRstd','FDmean','GSmean','GSstd','RVmean','RVstd','PPGmean','PPGstd')

%%  ---------------------------------------------------

X = [];
x = HRmean;  x = reshape(x',[4 100])'; x = mean(x,2);  X = [X,x];
x = HRstd;  x = reshape(x',[4 100])'; x = mean(x,2);  X = [X,x];
x = FDmean;  x = reshape(x',[4 100])'; x = mean(x,2);  X = [X,x];
x = GSmean;  x = reshape(x',[4 100])'; x = mean(x,2);  X = [X,x];
x = GSstd;  x = reshape(x',[4 100])'; x = mean(x,2);  X = [X,x];
x = RVmean;  x = reshape(x',[4 100])'; x = mean(x,2);  X = [X,x];
x = RVstd;  x = reshape(x',[4 100])'; x = mean(x,2);  X = [X,x];
x = PPGmean;  x = reshape(x',[4 100])'; x = mean(x,2);  X = [X,x];
x = PPGstd;  x = reshape(x',[4 100])'; x = mean(x,2);  X = [X,x];

FC = corr(X);  figure, imagesc(FC)

    

%%    fMRI feature for eacj subject

x = r_all(:,3);
xSubj = reshape(x',[4 100])';

x = squeeze(FC_partial_all(:,  3  ,4));
% x = squeeze(FC_all(:,  3  ,4));
xSubj = reshape(x',[4 100])';

ICC(xSubj,'C-1')
xfMRI = mean(xSubj,2);


%%  --------------------------------------------------

% xfMRI = X(:,5);
indConf = [2]  ;
[~,nCol] = size(BehMeas_num);


r_Beh = zeros(nCol,1);
p_Beh = zeros(nCol,1);
for c = 1:nCol
    Beh = BehMeas_num{:,c};    
    
    [r_tmp p_tmp]  = corr(Beh,xfMRI,'rows','complete');    
%     [r_tmp p_tmp]  = partialcorr(Beh,xfMRI,X(:,indConf),'rows','complete')

    
    r_Beh(c) = r_tmp;
    p_Beh(c) = p_tmp;
end

% plot(r_Beh)

ind_signif =  find(p_Beh<0.05);

BehMeas_categ = BehMeas_num.Properties.VariableNames;

    figure('position',[1000         341        1223         997])
for i = 1:length(ind_signif)
    subplot(4,3,i)
    c= ind_signif(i);
    Beh = BehMeas_num{:,c};
    
    [r_tmp p_tmp]  = corr(Beh,xfMRI,'rows','complete');
%         [r_tmp p_tmp]  = partialcorr(Beh,xfMRI,X(:,indConf),'rows','complete')
    
    scatter(Beh, xfMRI)
    h = refline;
    h.Color = 'r';
    ylabel('Memory of DCPM')
    xlabel(BehMeas_categ{c},'Interpreter','none')
    title(sprintf('r=%3.2f; p=%3.4f ', r_tmp, p_tmp))

end


%%     Hematocrit vs r_PPG

% tmp = r_all(:,4);
tmp = squeeze(FC_partial_all(:,  3  ,4));
tmpSubj = reshape(tmp',[4 100])';
xfMRI = mean(tmpSubj,2);

xfMRI = X(:,8);

%         [r_tmp p_tmp]  = partialcorr(Beh,xfMRI,X(:,indConf),'rows','complete')

 y = xfMRI; 
 
x = hematocrit;
% x = PressurePulse;
% y = X(:,8);
% x = X(:,7)

scatter(x,y,50,'k','LineWidth',2)
[r p]  = corr(x,y,'rows','complete')

model = fitlm(x,y)
pVal = model.Coefficients.pValue(2);
ax=gca; ax.FontWeight='bold'; hold on
Par = model.Coefficients.Estimate; P12 = Par(1)+Par(2)*1;
x_scal = (min(x)-0.2:0.0001:max(x)+0.2);
y_scal = Par(1) + Par(2)*x_scal;
plot(x_scal,y_scal,'r','LineWidth',3)

% xlabel('Hematocrit (%)')
xlabel('Heart rate variability (bpm)')
xlabel('Standard deviation of RV (bpm)')
ylabel('Partial correlation (GS,X_{RV})')
grid on






