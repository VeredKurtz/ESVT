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

%%
color_pareto = [0 112 192]/255;
color_uniform = [191 32 37]/255;
color_pareto_scatter = [147 209 255]/255;
color_uniform_scatter = [238 156 158]/255;

%% sub 1603 (risk neutral)
sub_rows_risk_neutral = (data(:,5)==1603);
data_risk_neutral = data(sub_rows_risk_neutral,:);

% UNIFORM dataset
uniform_rows_risk_neutral = find(data_risk_neutral(:,45)==0);
sv_uniform(:,1) = data_risk_neutral(uniform_rows_risk_neutral,20);
sv_uniform(:,2) = data_risk_neutral(uniform_rows_risk_neutral,22);

x1_high_uniform =  data_risk_neutral(uniform_rows_risk_neutral,8);
x1_low_uniform = data_risk_neutral(uniform_rows_risk_neutral,9);

% 2d hist
figure;
hist3(sv_uniform,[25 25],'CdataMode','auto')
xlabel('v1'); ylabel('v2');
colorbar; view(2); clim([0 4]);
% exportgraphics(gcf,'uniform2d_RN.eps','BackgroundColor','none','ContentType','vector');

% 1d hist
figure;
histogram(sv_uniform(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]);
% exportgraphics(gcf,'sv2_uniform_dist_RN.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - high outcome
figure;
histogram(x1_high_uniform(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[236 125 48]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_high_uniform_dist_RN.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - low outcome
figure;
histogram(x1_low_uniform(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[24 112 185]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_low_uniform_dist_RN.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - EV
ev_risk_neutral_uniform = 0.5.*x1_high_uniform + 0.5.*x1_low_uniform; 
figure;
histogram(ev_risk_neutral_uniform(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[0.7 0.7 0.7], 'FaceAlpha',1); ylim([0 0.3]);
exportgraphics(gcf,'ev_uniform_dist_RN.eps','BackgroundColor','none','ContentType','vector');

%%
% PARETO dataset
pareto_rows_risk_neutral = find(data_risk_neutral(:,45)==1);
sv_pareto(:,1) = data_risk_neutral(pareto_rows_risk_neutral,20);
sv_pareto(:,2) = data_risk_neutral(pareto_rows_risk_neutral,22);

x1_high_pareto =  data_risk_neutral(pareto_rows_risk_neutral,8);
x1_low_pareto = data_risk_neutral(pareto_rows_risk_neutral,9);

% 2d hist
figure;
hist3(sv_pareto,[25 25],'CdataMode','auto')
xlabel('v1'); ylabel('v2');
colorbar; view(2); clim([0 4]);
% exportgraphics(gcf,'pareto2d_RN.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - SV
figure;
histogram(sv_pareto(:,2),[0:4.5:60],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]);
exportgraphics(gcf,'sv2_pareto_dist_RN.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - high outcome
figure;
histogram(x1_high_pareto(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[236 125 48]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_high_pareto_dist_RN.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - low outcome
figure;
histogram(x1_low_pareto(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[24 112 185]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_low_pareto_dist_RN.eps','BackgroundColor','none','ContentType','vector');


% 1d hist - EV
ev_risk_neutral_pareto = 0.5.*x1_high_pareto + 0.5.*x1_low_pareto; 
figure;
histogram(ev_risk_neutral_pareto(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]);
exportgraphics(gcf,'ev_pareto_dist_RN.eps','BackgroundColor','none','ContentType','vector');


%% sub 1011 (risk seeking)
sub_rows_risk_seeking = (data(:,5)==1011);
data_risk_seeking = data(sub_rows_risk_seeking,:);

% UNIFORM dataset
uniform_rows_risk_seeking = find(data_risk_seeking(:,45)==0);
sv_uniform(:,1) = data_risk_seeking(uniform_rows_risk_seeking,20);
sv_uniform(:,2) = data_risk_seeking(uniform_rows_risk_seeking,22);

x1_high_uniform =  data_risk_seeking(uniform_rows_risk_seeking,8);
x1_low_uniform = data_risk_seeking(uniform_rows_risk_seeking,9);

% 2d hist
figure;
hist3(sv_uniform,[25 25],'CdataMode','auto')
xlabel('v1'); ylabel('v2');
colorbar; view(2); clim([0 4]);
% exportgraphics(gcf,'uniform2d_RS.eps','BackgroundColor','none','ContentType','vector');

% 1d hist
figure;
histogram(sv_uniform(:,1),[0:55:732],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]);
% exportgraphics(gcf,'sv2_uniform_dist_RS.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - high outcome
figure;
histogram(x1_high_uniform(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[236 125 48]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_high_uniform_dist_RS.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - low outcome
figure;
histogram(x1_low_uniform(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[24 112 185]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_low_uniform_dist_RS.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - EV
ev_risk_seeking_uniform = 0.5.*x1_high_uniform + 0.5.*x1_low_uniform; 
figure;
histogram(ev_risk_seeking_uniform(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]);
exportgraphics(gcf,'ev_uniform_dist_RS.eps','BackgroundColor','none','ContentType','vector');

%% PARETO dataset
pareto_rows_risk_seeking = find(data_risk_seeking(:,45)==1);
sv_pareto(:,1) = data_risk_seeking(pareto_rows_risk_seeking,20);
sv_pareto(:,2) = data_risk_seeking(pareto_rows_risk_seeking,22);

x1_high_pareto =  data_risk_seeking(pareto_rows_risk_seeking,8);
x1_low_pareto = data_risk_seeking(pareto_rows_risk_seeking,9);

% 2d hist
figure;
hist3(sv_pareto,[25 25],'CdataMode','auto')
xlabel('v1'); ylabel('v2');
colorbar; view(2); clim([0 4]);
% exportgraphics(gcf,'pareto2d_RS.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - SV
figure;
histogram(sv_pareto(:,2),[0:55:730],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]);
% exportgraphics(gcf,'sv2_pareto_dist_RS.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - high outcome
figure;
histogram(x1_high_pareto(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[236 125 48]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_high_pareto_dist_RS.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - low outcome
figure;
histogram(x1_low_pareto(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[24 112 185]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_low_pareto_dist_RS.eps','BackgroundColor','none','ContentType','vector');


% 1d hist - EV
ev_risk_seeking_pareto = 0.5.*x1_high_pareto + 0.5.*x1_low_pareto; 
figure;
histogram(ev_risk_seeking_pareto(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]);
exportgraphics(gcf,'ev_pareto_dist_RS.eps','BackgroundColor','none','ContentType','vector');


%% sub 1619 (risk averse)
sub_rows_risk_averse = (data(:,5)==1619);
data_risk_averse = data(sub_rows_risk_averse,:);

% UNIFORM dataset
uniform_rows_risk_averse = find(data_risk_averse(:,45)==0);
sv_uniform(:,1) = data_risk_averse(uniform_rows_risk_averse,20);
sv_uniform(:,2) = data_risk_averse(uniform_rows_risk_averse,22);

x1_high_uniform =  data_risk_averse(uniform_rows_risk_averse,8);
x1_low_uniform = data_risk_averse(uniform_rows_risk_averse,9);

% 2d hist
figure;
hist3(sv_uniform,[25 25],'CdataMode','auto')
xlabel('v1'); ylabel('v2');
colorbar; view(2); clim([0 4]);
% exportgraphics(gcf,'uniform2d_RA.eps','BackgroundColor','none','ContentType','vector');

% 1d hist
figure;
histogram(sv_uniform(:,1),[1:0.4:6.5],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]);
% exportgraphics(gcf,'sv2_uniform_dist_RA.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - high outcome
figure;
histogram(x1_high_uniform(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[236 125 48]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_high_uniform_dist_RA.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - low outcome
figure;
histogram(x1_low_uniform(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[24 112 185]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_low_uniform_dist_RA.eps','BackgroundColor','none','ContentType','vector');


% 1d hist - EV
ev_risk_averse_uniform = 0.5.*x1_high_uniform + 0.5.*x1_low_uniform; 
figure;
histogram(ev_risk_averse_uniform(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]);
exportgraphics(gcf,'ev_uniform_dist_RA.eps','BackgroundColor','none','ContentType','vector');

%% PARETO dataset
pareto_rows_risk_averse = find(data_risk_averse(:,45)==1);
sv_pareto(:,1) = data_risk_averse(pareto_rows_risk_averse,20);
sv_pareto(:,2) = data_risk_averse(pareto_rows_risk_averse,22);

x1_high_pareto =  data_risk_averse(pareto_rows_risk_averse,8);
x1_low_pareto = data_risk_averse(pareto_rows_risk_averse,9);

% 2d hist
figure;
hist3(sv_pareto,[25 25],'CdataMode','auto')
xlabel('v1'); ylabel('v2');
colorbar; view(2); clim([0 4]);
% exportgraphics(gcf,'pareto2d_RA.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - SV
figure;
histogram(sv_pareto(:,2),[0:0.45:6],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]); xlim([0.2 6]);
% exportgraphics(gcf,'sv2_pareto_dist_RA.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - high outcome
figure;
histogram(x1_high_pareto(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[236 125 48]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_high_pareto_dist_RA.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - low outcome
figure;
histogram(x1_low_pareto(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[24 112 185]./255, 'FaceAlpha',1); ylim([0 0.3]);
% exportgraphics(gcf,'x1_low_pareto_dist_RA.eps','BackgroundColor','none','ContentType','vector');

% 1d hist - EV
ev_risk_averse_pareto = 0.5.*x1_high_pareto + 0.5.*x1_low_pareto; 
figure;
histogram(ev_risk_averse_pareto(:,1),[0:4.5:60],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); ylim([0 0.3]);
exportgraphics(gcf,'ev_pareto_dist_RA.eps','BackgroundColor','none','ContentType','vector');
