

load_path_RETR = ['Export/2019_11_01/CP_Model_order_selection_RETR/']; 
load([load_path_RETR,'RETROICOR_random_voxels'])
r_RETR_all = r_all;
r_RETR_all = squeeze(mean(r_RETR_all,3));

load_path_DCPM = ['Export/2019_11_01/CP_Model_order_selection_DCPM/']; 
load([load_path_DCPM,'CPM_random_voxels'])
r_DCPM_all = r_all;
r_DCPM_all = squeeze(mean(r_DCPM_all,3));

% Better upload workspace 

%%

[nScans, nLags, nVoxels, nModels] = size(r_all);  nSubj = nScans/4;

r_RETR_sc = squeeze(mean(r_RETR_all,1));
r_DCPM_sc = squeeze(mean(r_DCPM_all,1));

plot(temp_shift_scale,r_DCPM_sc), grid on
xlabel('Time lag (s)')
ylabel('Mean correlation')

for i = 1:nModels
    legend_labels{i} = sprintf('M.O.: %d', i)
end
legend(legend_labels)
% title('Dynamic cardiac pulsatility model (DCPM)')
title('RETROICOR')


%%    Average first within subject

model_1 = squeeze(r_RETR_all(:,25,6));
model_2 = squeeze(r_DCPM_all(:,30,6));

model_1_subj = zeros(nSubj,1);
model_2_subj = zeros(nSubj,1);
for i = 1 :nSubj
    ind = (1:4) + (i-1)*4;
    model_1_subj(i) = mean(model_1(ind));
    model_2_subj(i) = mean(model_2(ind));
end

fprintf('Mean correlation: %3.3f+-%3.3f     \n',mean(model_1_subj), std(model_1_subj))
fprintf('Mean correlation: %3.3f+-%3.3f     \n',mean(model_2_subj), std(model_2_subj))

[ttest_p ttest_h] =  ttest(model_1_subj, model_2_subj)




%% Improvement with higher order  when at optimal time lag   ----------------

r_RETR_t_opt= squeeze(r_RETR_all(:,25,:));
r_DCPM_t_opt = squeeze(r_DCPM_all(:,30,:));


model_order = 1:8;
r_mean_RETR = mean(r_RETR_t_opt); r_std_RETR = std(r_RETR_t_opt)/sqrt(nSubj);
r_mean_DCPM = mean(r_DCPM_t_opt); r_std_DCPM = std(r_DCPM_t_opt)/sqrt(nSubj);

figure
shadedErrorBar(model_order,r_mean_RETR,r_std_RETR, 'ro') , hold on
shadedErrorBar(model_order,r_mean_DCPM,r_std_DCPM, 'ko',0.9) , hold on

grid on, xlim([1 nModels])
xlabel('Model order')
ylabel('Mean Correlation')
title('Mean correlation at optimal time lag')







