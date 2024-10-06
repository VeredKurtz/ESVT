clear
clc

file_name = 'ESVT_data_allSubjects.xlsx';
path = '/Users/kurtzv02/Dropbox/DN - uniform,pareto,2,6/data analysis codes/Vered/';
[data, titles, raw] = xlsread([path file_name]);

% drop subs 1511 and 1605, who are missing the cool down trials
sub1511_rows = find(data(:,5)==1511);
data(sub1511_rows,:) = [];

sub1605_rows = find(data(:,5)==1605);
data(sub1605_rows,:) = [];

%%
color_pareto = [0 112 192]/255;
color_uniform = [191 32 37]/255;
color_pareto_scatter = [147 209 255]/255;
color_uniform_scatter = [238 156 158]/255;


%% drop 6-options
num_options = data(:,1);
six_options_rows = find(num_options==6);
data(six_options_rows,:) = [];

%% important parameters
riskAlpha = data(:,6);
ChoiceMaxSV = data(:,37); 
Pareto_dummy = data(:,45);
SV_option1 = data(:,20);
SV_option2 = data(:,22);
cool_down_dummy = data(:,35);

%%
sid = data(:,5);
subjects = unique(sid);
maxSV = data(:,36);

% prepare a matrix for the results
riskAlpha_sllSubs = NaN(length(subjects),1);

for sub = 1:length(subjects)
    subNum = subjects(sub);
    sub_rows = find(sid==subNum);
    sub_riskAlpha = riskAlpha(sub_rows);
    riskAlpha_sllSubs(sub,1) = sub_riskAlpha(1,1);
    clear sub_rows sub_riskAlpha subNum
end

%%
figure;
X=0:1:4;
histogram(riskAlpha_sllSubs,15,'Normalization','probability','FaceColor',[0.5 0.5 0.5],'LineWidth',1.5);
ylim([0 0.25]); xlabel('risk preferences parameters (\rho)','FontSize',12); ylabel('frequency (%)','FontSize',12); set(gca,'FontSize',12);
exportgraphics(gcf,'distRho.eps','BackgroundColor','none','ContentType','vector');
