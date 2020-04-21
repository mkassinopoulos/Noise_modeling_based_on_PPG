


path_load = 'E:\CloudStation\HCP\Analysis\Pulsatility_model\Export\2020_02_01\CP_Model_order_selection\';

flag_signiciant = 0 ;

if flag_signiciant == 1
    load([path_load,'r_all_Cardiac_RETROICOR.mat']); r_RETR = r_all;
    load([path_load,'r_all_Cardiac_DCPM_CA.mat']); r_CPM_CA = r_all;
    load([path_load,'r_all_Cardiac_DCPM_VA.mat']); r_CPM_VA = r_all;
else
    load([path_load,'r_all_random_RETROICOR.mat']); r_RETR = r_all;
    load([path_load,'r_all_random_DCPM_CA.mat']); r_CPM_CA = r_all;
    load([path_load,'r_all_random_DCPM_VA.mat']); r_CPM_VA = r_all;
end

%%  -------------------------

r_all = r_RETR;
tmp = squeeze(mean(r_all,3));
r_tmp_v1 = squeeze(mean(tmp,1));
max(r_tmp_v1(:))

figure
c = parula(8); % c = brewermap(8,'PuRd');
for i  = 1:8
    plot(-temp_shift_scale,r_tmp_v1(:,i),'Color',c(i,:)), hold on
end

grid, xlabel('Lag time (s)'), ylabel('Mean correlation')

title('RETROICOR')
title('CPM_{CA}')
title('CPM_{VA}')

ylim([0.14  0.29])
ylim([0.04  0.1])


%%  ----------------------------

model_order = 1:8;
nSubj = 100;

ind_lag =  25; 
order = squeeze(r_RETR(:,ind_lag,:,:));
order = squeeze(mean(order,2));
order_sbj = zeros(nSubj,8);
for i = 1:nSubj
    ind = (i-1)*4 + (1:4);
    x = order(ind,:);
    order_sbj(i,:) = mean(x);
end
r_mean_RETR = mean(order_sbj);
r_std_RETR = std(order_sbj)/sqrt(nSubj);


ind_lag =  30; 
order = squeeze(r_CPM_CA(:,ind_lag,:,:));
order = squeeze(mean(order,2));
order_sbj = zeros(nSubj,8);
for i = 1:nSubj
    ind = (i-1)*4 + (1:4);
    x = order(ind,:);
    order_sbj(i,:) = mean(x);
end
RETR_order = order_sbj;
r_mean_CPM_CA = mean(order_sbj);
r_std_CPM_CA = std(order_sbj)/sqrt(nSubj);


ind_lag =  30; 
order = squeeze(r_CPM_VA(:,ind_lag,:,:));
order = squeeze(mean(order,2));
order_sbj = zeros(nSubj,8);
for i = 1:nSubj
    ind = (i-1)*4 + (1:4);
    x = order(ind,:);
    order_sbj(i,:) = mean(x);
end
RETR_order = order_sbj;
r_mean_CPM_VA = mean(order_sbj);
r_std_CPM_VA = std(order_sbj)/sqrt(nSubj);


shadedErrorBar(model_order,r_mean_RETR,r_std_RETR, 'ro') , hold on
shadedErrorBar(model_order,r_mean_CPM_CA,r_std_CPM_CA, 'ko',0.9) , hold on
shadedErrorBar(model_order,r_mean_CPM_VA,r_std_CPM_VA, 'go',0.9) , hold on


grid on, xlim([1 8.1])
xlabel('Model order')
ylabel('Mean correlation')
title('Mean correlation at optimal lag time')

% ylim([0.19 0.30])
ylim([0.065 0.105])


%%  Compare models - p-value for optimal order and shift  -----------------------

%  RETR(M=2, dt=0) --> RETR(M=2, dt=0.4)

ind_lag =  21; 
r_opt = squeeze(r_RETR(:,ind_lag,:,2));
r_opt = mean(r_opt,2);
r_opt = reshape(r_opt',[4 100])';
r_opt_x1 = mean(r_opt,2);
mean(r_opt_x1) , std(r_opt_x1)

ind_lag =  25; 
r_opt = squeeze(r_RETR(:,ind_lag,:,2));
r_opt = mean(r_opt,2);
r_opt = reshape(r_opt',[4 100])';
r_opt_x2 = mean(r_opt,2);
mean(r_opt_x2)

[ttest_p ttest_h] = ttest(r_opt_x1, r_opt_x2)


%  RETR(M=2, dt=0.4) --> RETR(M=6, dt=0.4)

ind_lag =  25; 
r_opt = squeeze(r_RETR(:,ind_lag,:,2));
r_opt = mean(r_opt,2);
r_opt = reshape(r_opt',[4 100])';
r_opt_x1 = mean(r_opt,2);
mean(r_opt_x1) , std(r_opt_x1)

ind_lag =  25; 
r_opt = squeeze(r_RETR(:,ind_lag,:,6));
r_opt = mean(r_opt,2);
r_opt = reshape(r_opt',[4 100])';
r_opt_x2 = mean(r_opt,2);
mean(r_opt_x2)

[ttest_p ttest_h] = ttest(r_opt_x1, r_opt_x2)



%  RETR(M=6, dt=0.4) --> CPM_CA(M=6, dt=0.9)

ind_lag =  25; 
r_opt = squeeze(r_RETR(:,ind_lag,:,6));
r_opt = mean(r_opt,2);
r_opt = reshape(r_opt',[4 100])';
r_opt_x1 = mean(r_opt,2);
mean(r_opt_x1) , std(r_opt_x1)

ind_lag =  30; 
r_opt = squeeze(r_CPM_CA(:,ind_lag,:,6));
r_opt = mean(r_opt,2);
r_opt = reshape(r_opt',[4 100])';
r_opt_x2 = mean(r_opt,2);
mean(r_opt_x2) , std(r_opt_x2)


[ttest_p ttest_h] = ttest(r_opt_x1, r_opt_x2)


%%   HRV vs dr


x1 =  squeeze(r_RETR(:,25,:,6)); x1 = mean(x1,2);
tmp = reshape(x1,[4 100])'; x1subj = mean(tmp,2);

x3 =  squeeze(r_CPM_CA(:,30,:,6)); x3 = mean(x3,2);
tmp = reshape(x3,[4 100])'; x3subj = mean(tmp,2);

[ttest_h, ttest_p] = ttest(x1subj,x3subj)

dx_sbj = x3subj-x1subj;
dx_relat = (x3subj-x1subj)*100./x1subj;

HRstd_sbj = reshape(HRstd_all,[4 100])'; HRstd_sbj = mean(HRstd_sbj,2);


%%  ------------------

 y = dx_relat;
%  y = dx_sbj
 
x = HRstd_sbj;
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

ylabel('% improvement')
% ylabel('?r(RETR-->CPM_{CA})')

grid on











