clear
clc

file_name = 'ESVT_data_allSubjects.xlsx';
path = '/Users/kurtzv02/Dropbox/DN - uniform,pareto,2,6/data analysis codes/Vered/';
[data, titles, raw] = xlsread([path file_name]);

% remove 6-options
six_options_rows = find(data(:,1)==6);
data(six_options_rows,:) = [];

% remove cool-down trials
coold_down_rows = find(data(:,35)==1);
data(coold_down_rows,:) = [];

% drop subject 1511 who did not complete both datasets
rows_1511 = find(data(:,5)==1511);
data(rows_1511,:) = [];

% drop subject 1605 who did not complete both datasets
rows_1605 = find(data(:,5)==1605);
data(rows_1605,:) = [];

% drop subject 911 who did not complete both datasets
rows_911 = find(data(:,5)==911);
data(rows_911,:) = [];

% drop subject 913 who did not complete both datasets
rows_913 = find(data(:,5)==913);
data(rows_913,:) = [];


%% colors
color_pareto = [0 112 192]/255;
color_uniform = [191 32 37]/255;
color_pareto_scatter = [147 209 255]/255;
color_uniform_scatter = [238 156 158]/255;


%% split by treatment

data_pareto_rows = find(data(:,45)==1);
data_pareto = data(data_pareto_rows,:);
data_uniform_rows = find(data(:,45)==0);
data_uniform = data(data_uniform_rows,:);

%%
sid = data(:,5);
subjects = unique(sid);
num_subjects = length(subjects);

%% subjects loop 

for sub=1:num_subjects
 current_sub = subjects(sub);

 % pareto
 sub_rows_pareto = find(data_pareto(:,5)==current_sub);
 ev_option_1_pareto = data_pareto(sub_rows_pareto,21);
 median_option1_pareto(sub,1) = nanmedian(ev_option_1_pareto);
 ev_option_2_pareto = data_pareto(sub_rows_pareto,23);
 sv_option_1_pareto = data_pareto(sub_rows_pareto,20);
 median_sv_option1_pareto(sub,1) = nanmedian(sv_option_1_pareto);
 sv_option_2_pareto = data_pareto(sub_rows_pareto,22);
 chosen_option_pareto = data_pareto(sub_rows_pareto,3);
 for trial = 1:length(sub_rows_pareto)
      if ev_option_1_pareto(trial)>=ev_option_2_pareto(trial)
         max_ev_pareto(trial) = 1;
     else
         max_ev_pareto(trial) = 2;
      end
      if sv_option_1_pareto(trial)>=sv_option_2_pareto(trial)
         max_sv_pareto(trial) = 1;
     else
         max_sv_pareto(trial) = 2;
     end
     if chosen_option_pareto(trial)==max_ev_pareto(trial)
         chose_max_ev_pareto(trial,1) = 1;
     else
         chose_max_ev_pareto(trial,1) = 0;
     end
     if chosen_option_pareto(trial)==max_sv_pareto(trial)
         chose_max_sv_pareto(trial,1) = 1;
     else
         chose_max_sv_pareto(trial,1) = 0;
     end

     if chosen_option_pareto(trial) == 1
         chosen_option_high_pareto = data_pareto(sub_rows_pareto(trial),8);
         chosen_option_low_pareto = data_pareto(sub_rows_pareto(trial),9);
         unchosen_option_high_pareto = data_pareto(sub_rows_pareto(trial),10);
         unchosen_option_low_pareto = data_pareto(sub_rows_pareto(trial),11);
         if chosen_option_high_pareto<unchosen_option_high_pareto && chosen_option_low_pareto<unchosen_option_low_pareto
             violateFOSD_pareto(trial,1) = 1;
         else
             violateFOSD_pareto(trial,1) = 0;
         end
     elseif chosen_option_pareto(trial) == 2
         chosen_option_high_pareto = data_pareto(sub_rows_pareto(trial),10);
         chosen_option_low_pareto = data_pareto(sub_rows_pareto(trial),11);
         unchosen_option_high_pareto = data_pareto(sub_rows_pareto(trial),8);
         unchosen_option_low_pareto = data_pareto(sub_rows_pareto(trial),9);
         if chosen_option_high_pareto<unchosen_option_high_pareto && chosen_option_low_pareto<unchosen_option_low_pareto
             violateFOSD_pareto(trial,1) = 1;
         else
             violateFOSD_pareto(trial,1) = 0;
         end
     end

 end
 results(sub,1) = mean(chose_max_ev_pareto);
 results(sub,2) = length(find(violateFOSD_pareto==1));
 results(sub,3) = mean(chose_max_sv_pareto);

 % uniform
 sub_rows_uniform = find(data_uniform(:,5)==current_sub);
 ev_option_1_uniform = data_uniform(sub_rows_uniform,21);
 median_option1_uniform(sub,1) = nanmedian(ev_option_1_uniform);
 ev_option_2_uniform = data_uniform(sub_rows_uniform,23);
 sv_option_1_uniform = data_uniform(sub_rows_uniform,20);
 median_sv_option1_uniform(sub,1) = nanmedian(sv_option_1_uniform);
 sv_option_2_uniform = data_uniform(sub_rows_uniform,22);
 chosen_option_uniform = data_uniform(sub_rows_uniform,3);

 %  stats
 diff_medians_sv = (median_sv_option1_uniform-median_sv_option1_pareto); 
 mean_median_uniform = nanmedian(median_sv_option1_uniform);
 mean_median_pareto = nanmedian(median_sv_option1_pareto);
 [p,h,stats] = signrank(median_sv_option1_uniform,median_sv_option1_pareto,"tail","right");

 for trial = 1:length(sub_rows_uniform)
     if ev_option_1_uniform(trial)>=ev_option_2_uniform(trial)
         max_ev_uniform(trial) = 1;
     else
         max_ev_uniform(trial) = 2;
     end
     if sv_option_1_uniform(trial)>=sv_option_2_uniform(trial)
         max_sv_uniform(trial) = 1;
     else
         max_sv_uniform(trial) = 2;
     end
     if chosen_option_uniform(trial)==max_ev_uniform(trial)
         chose_max_ev_uniform(trial,1) = 1;
     else
         chose_max_ev_uniform(trial,1) = 0;
     end
     if chosen_option_uniform(trial)==max_sv_uniform(trial)
         chose_max_sv_uniform(trial,1) = 1;
     else
         chose_max_sv_uniform(trial,1) = 0;
     end

     if chosen_option_uniform(trial) == 1
         chosen_option_high_uniform = data_uniform(sub_rows_uniform(trial),8);
         chosen_option_low_uniform = data_uniform(sub_rows_uniform(trial),9);
         unchosen_option_high_uniform = data_uniform(sub_rows_uniform(trial),10);
         unchosen_option_low_uniform = data_uniform(sub_rows_uniform(trial),11);
         if chosen_option_high_uniform<unchosen_option_high_uniform && chosen_option_low_uniform<unchosen_option_low_uniform
             violateFOSD_uniform(trial,1) = 1;
         else
             violateFOSD_uniform(trial,1) = 0;
         end
     elseif chosen_option_uniform(trial) == 2
         chosen_option_high_uniform = data_uniform(sub_rows_uniform(trial),10);
         chosen_option_low_uniform = data_uniform(sub_rows_uniform(trial),11);
         unchosen_option_high_uniform = data_uniform(sub_rows_uniform(trial),8);
         unchosen_option_low_uniform = data_uniform(sub_rows_uniform(trial),9);
         if chosen_option_high_uniform<unchosen_option_high_uniform && chosen_option_low_uniform<unchosen_option_low_uniform
             violateFOSD_uniform(trial,1) = 1;
         else
             violateFOSD_uniform(trial,1) = 0;
         end
     end

 end
 results(sub,4) = mean(chose_max_ev_uniform);
 results(sub,5) = length(find(violateFOSD_uniform==1));
 results(sub,6) = mean(chose_max_sv_uniform);

 clear chose_max_ev_pareto chosen_option_pareto max_ev_pareto ev_option_1_pareto ev_option_2_pareto sub_rows_pareto violateFOSD_pareto
 clear chose_max_ev_uniform chosen_option_uniform max_ev_uniform ev_option_1_uniform ev_option_2_uniform sub_rows_uniform violateFOSD_uniform

end

%% figures
zero_locations_pareto = find(results(:,2)==0);
zero_locations_uniform = find(results(:,5)==0);
results(zero_locations_pareto,2) = 0.0001;
results(zero_locations_uniform,5) = 0.0001;

figure;
subplot(1,2,1);
maxEV_Table = array2table(100.*[results(:,1) results(:,4)],...
                'VariableNames',{'Pareto','Uniform'});

v1=violinplot(maxEV_Table, [] ,'Width', 0.3,'ViolinColor',[color_pareto; color_uniform], 'Bandwidth',5);
ylim([0 100]); xlim([0.5 2.5]); set(gca,'FontSize',14); xticks([1 2]); xticklabels({'Pareto','Uniform'});
ylabel('chose max. EV (% of trials)');

subplot(1,2,2);
violateFOSD_Table = array2table([results(:,2) results(:,5)],...
                'VariableNames',{'Pareto','Uniform'});

v1=violinplot(violateFOSD_Table, [] ,'Width', 0.3,'ViolinColor',[color_pareto; color_uniform], 'Bandwidth',5);
ylim([0 150]); xlim([0.5 2.5]); set(gca,'FontSize',14); xticks([1 2]); xticklabels({'Pareto','Uniform'});
ylabel('num. of FOSD violations (#)');

figure;
subplot(1,2,1);
maxEV_Table = array2table(100.*[results(:,3) results(:,6)],...
                'VariableNames',{'Pareto','Uniform'});

v1=violinplot(maxEV_Table, [] ,'Width', 0.3,'ViolinColor',[color_pareto; color_uniform], 'Bandwidth',5);
ylim([0 100]); xlim([0.5 2.5]); set(gca,'FontSize',14); xticks([1 2]); xticklabels({'Pareto','Uniform'});
ylabel('chose max. SV (% of trials)');

subplot(1,2,2);
violateFOSD_Table = array2table([results(:,2) results(:,5)],...
                'VariableNames',{'Pareto','Uniform'});

v1=violinplot(violateFOSD_Table, [] ,'Width', 0.3,'ViolinColor',[color_pareto; color_uniform], 'Bandwidth',5);
ylim([0 150]); xlim([0.5 2.5]); set(gca,'FontSize',14); xticks([1 2]); xticklabels({'Pareto','Uniform'});
ylabel('num. of FOSD violations (#)');



