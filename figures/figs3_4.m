clear
clc

file_name = 'subject_level_fits_incl_noise_05102024.xlsx';
path = '/Users/kurtzv02/Dropbox/DN - uniform,pareto,2,6/results/';
[data, titles, raw] = xlsread([path file_name]);

%%
color_pareto = [0 112 192]/255;
color_uniform = [191 32 37]/255;
color_pareto_scatter = [147 209 255]/255;
color_uniform_scatter = [238 156 158]/255;

%% fig 1 - r (RUM) vs risk preference rho
subjects_rho = data(:,2);
subjects_r_uniform = data(:,13);
subjects_r_pareto = data(:,17);

% uniform
figure;
scatter(subjects_rho,subjects_r_uniform,160,'MarkerFaceColor',color_uniform_scatter,'MarkerEdgeColor','k','LineWidth',0.5);
hold on
x= 0:0.001:2.5;
y= ones(length(x),1);
y2 = 1./(x);
plot(x,y2,'LineWidth',1.5,'LineStyle','--','Color',color_uniform); ylim([0,3.5]); xlabel('\rho (STAGE I)','FontSize',12); ylabel('r parameter (DN, STAGE II)','FontSize',12);
opt_rho_curvature = 1./subjects_rho;
rmse1 = rmse(opt_rho_curvature,subjects_r_uniform);
 % exportgraphics(gcf,'rho_vs_r_uniform.eps','BackgroundColor','none','ContentType','vector');

% pareto
figure;
scatter(subjects_rho,subjects_r_pareto,160,'MarkerFaceColor',color_pareto_scatter,'MarkerEdgeColor','k','LineWidth',0.5);
hold on
x= 0:0.001:2.5;
y= ones(length(x),1);
y2 = 1./(x);
plot(x,y2,'LineWidth',1.5,'LineStyle','--','Color',color_pareto); ylim([0,3.5]); xlabel('\rho (STAGE I)','FontSize',12); ylabel('r parameter (DN, STAGE II)','FontSize',12); 
rmse2 = rmse(opt_rho_curvature,subjects_r_pareto);
 % exportgraphics(gcf,'rho_vs_r_pareto.eps','BackgroundColor','none','ContentType','vector');

%% alpha (ESVT) vs. risk preferences rho
figure;
alpha_uniform_ESVT = data(:,4);
scatter(subjects_rho,alpha_uniform_ESVT,160,'MarkerFaceColor',color_uniform_scatter,'MarkerEdgeColor','k','LineWidth',0.5);
hold on
x= 0:0.001:2.5;
y= 1./x;
plot(x,y,'LineWidth',1.5,'LineStyle','--','Color',color_uniform); ylim([0,6]); xlabel('\rho (STAGE I)','FontSize',12); ylabel('\alpha parameter (DN, STAGEII)','FontSize',12);
rmse3 = rmse(opt_rho_curvature,alpha_uniform_ESVT); 
% exportgraphics(gcf,'rho_vs_alpha_uniform.eps','BackgroundColor','none','ContentType','vector');

figure;
alpha_pareto_ESVT = data(:,9);
scatter(subjects_rho,alpha_pareto_ESVT,160,'MarkerFaceColor',color_pareto_scatter,'MarkerEdgeColor','k','LineWidth',0.5);
hold on
x= 0:0.001:2.5;
y= 1./(x);
plot(x,y,'LineWidth',1.5,'LineStyle','--','Color',color_pareto); ylim([0,6]); xlabel('\rho (STAGE I)','FontSize',12); ylabel('\alpha parameter (DN, STAGE II)','FontSize',12); 
rmse4 = rmse(opt_rho_curvature,alpha_pareto_ESVT); 
% exportgraphics(gcf,'rho_vs_alpha_pareto.eps','BackgroundColor','none','ContentType','vector');

%% fig. 2 - BIC RUM vs ESVT
% uniform 
BIC_RUM_uniform = data(:,16);
BIC_ESVT_uniform = data(:,7);
num_subjects = length(data(:,1));


figure;
scatter(BIC_RUM_uniform,BIC_ESVT_uniform,400,'MarkerFaceColor',color_uniform_scatter,'MarkerEdgeColor','k','LineWidth',0.5); 
hold on
x2 = 0:0.1:5000;
y2=x2;
plot(x2,y2,'LineWidth',1.5,'LineStyle','--','Color',color_uniform); 
xlim([0 750]); ylim([0 750]);
hold on
errorbar(mean(BIC_RUM_uniform),mean(BIC_ESVT_uniform),std(BIC_RUM_uniform)./sqrt(num_subjects),'horizontal','CapSize',0,'LineWidth',3,'Color','k');
hold on
errorbar(mean(BIC_RUM_uniform),mean(BIC_ESVT_uniform),std(BIC_ESVT_uniform)./sqrt(num_subjects),'vertical','CapSize',0,'LineWidth',3,'Color','k');
[p1,h1,stats1] = signrank(BIC_RUM_uniform,BIC_ESVT_uniform,'tail','right');
% exportgraphics(gcf,'BIC_RUM_ESVT_uniform.eps','BackgroundColor','none','ContentType','vector');

% outliers inset 
% figure;
% scatter(BIC_RUM_uniform,BIC_ESVT_uniform,160,'MarkerFaceColor',color_uniform_scatter,'MarkerEdgeColor','k','LineWidth',0.5); 
% hold on
% x2 = 0:0.1:5000;
% y2=x2;
% plot(x2,y2,'LineWidth',1.5,'LineStyle','--','Color',color_uniform); 
% xlim([1500 5000]); ylim([0 1500]);
% % exportgraphics(gcf,'outliers_BIC_uniform.eps','BackgroundColor','none','ContentType','vector');

 
% pareto
BIC_RUM_pareto = data(:,20);
BIC_ESVT_pareto = data(:,12);

figure;
scatter(BIC_RUM_pareto,BIC_ESVT_pareto,400,'MarkerFaceColor',color_pareto_scatter,'MarkerEdgeColor','k','LineWidth',0.5); 
hold on
x2 = 0:0.1:5000;
y2=x2;
plot(x2,y2,'LineWidth',1.5,'LineStyle','--','Color',color_pareto);
xlim([0 1500]); ylim([0 1500]);
hold on
errorbar(mean(BIC_RUM_pareto),mean(BIC_ESVT_pareto),std(BIC_RUM_pareto)./sqrt(num_subjects),'horizontal','CapSize',0,'LineWidth',3,'Color','k');
hold on
errorbar(mean(BIC_RUM_pareto),mean(BIC_ESVT_pareto),std(BIC_ESVT_pareto)./sqrt(num_subjects),'vertical','CapSize',0,'LineWidth',3,'Color','k');
% exportgraphics(gcf,'BIC_RUM_ESVT_pareto.eps','BackgroundColor','none','ContentType','vector');
[p2,h2,stats2] = signrank(BIC_RUM_pareto,BIC_ESVT_pareto,'tail','right');

% outliers inset 
% figure;
% scatter(BIC_RUM_pareto,BIC_ESVT_pareto,160,'MarkerFaceColor',color_pareto_scatter,'MarkerEdgeColor','k','LineWidth',0.5); 
% hold on
% x2 = 0:0.1:1500;
% y2=x2;
% plot(x2,y2,'LineWidth',1.5,'LineStyle','--','Color',color_pareto);
% xlim([750 1500]); ylim([0 1500]);
% exportgraphics(gcf,'outliers_BIC_pareto.eps','BackgroundColor','none','ContentType','vector');


%% insets - histograms of the differences
% uniform
BIC_diff_uniform = BIC_RUM_uniform-BIC_ESVT_uniform;
edges2 = [-600:50:600];
figure;
histogram(BIC_diff_uniform,edges2,'FaceColor',color_uniform_scatter,'EdgeColor',color_uniform); xlabel('\Delta BIC'); ylabel('num. of subjects'); ylim([0 40]); xlim([-600 600]);
% exportgraphics(gcf,'BIC_diff_uniform.eps','BackgroundColor','none','ContentType','vector');

% pareto
BIC_diff_pareto = BIC_RUM_pareto-BIC_ESVT_pareto;
figure;
histogram(BIC_diff_pareto,edges2,'FaceColor',color_pareto_scatter,'EdgeColor',color_pareto); xlabel('\Delta BIC'); ylabel('num. of subjects'); ylim([0 40]); xlim([-600 600]);
% exportgraphics(gcf,'BIC_diff_pareto.eps','BackgroundColor','none','ContentType','vector');

%% fig. 3 - M - uinform vs pareto (ESVT)
M_uniform = data(:,3);
M_pareto = data(:,8);

figure;
scatter(M_uniform,M_pareto,300,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k','LineWidth',0.5); 
 xlim([0 100]); ylim([0 100]); 
hold on
x3 = 0:0.01:100;
y3=x3;
plot(x3,y3,'LineWidth',1.5,'LineStyle','--','Color','K');
 hold on
 errorbar(mean(M_uniform(M_uniform<100)),mean(M_pareto(M_pareto<100)),std(M_uniform(M_uniform<100))./sqrt(num_subjects),'horizontal','CapSize',0,'LineWidth',3,'Color','k');
 hold on
 errorbar(mean(M_uniform(M_uniform<100)),mean(M_pareto(M_pareto<100)),std(M_pareto(M_pareto<100))./sqrt(num_subjects),'vertical','CapSize',0,'LineWidth',3,'Color','k');
 % exportgraphics(gcf,'M_diff.eps','BackgroundColor','none','ContentType','vector');

% outliers inset 
figure;
scatter(M_uniform,M_pareto,160,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k','LineWidth',0.5); 
xlim([100 400]); ylim([100 400]); 
hold on
x4 = 100:0.01:1000;
y4=x4;
plot(x4,y4,'LineWidth',1.5,'LineStyle','--','Color','K'); 
 % exportgraphics(gcf,'M_diff_inset.eps','BackgroundColor','none','ContentType','vector');

% histogram inset
figure;
M_diff = M_uniform - M_pareto;
histogram(M_diff(M_diff<200),20,'FaceColor',[0.7 0.7 0.7],'EdgeColor','k'); xlabel('\Delta M'); ylabel('num. of subjects'); ylim([0 45]); xlim([-200 200]);
% exportgraphics(gcf,'M_diff_histogram.eps','BackgroundColor','none','ContentType','vector');

[p3,h3,stats3] = signrank(M_uniform,M_pareto,'tail','right');

%% fig. 4 - r distributions (RUM)
% edges = 0:0.175:3.5;
% figure;
% histogram(subjects_r_pareto,edges,'FaceColor',color_pareto_scatter,'EdgeColor',color_pareto,'FaceAlpha',1,'Orientation','horizontal');
% xlim([0,20]);
% % exportgraphics(gcf,'r_dist_RUM_pareto.eps','BackgroundColor','none','ContentType','vector');
% 
% figure;
% histogram(subjects_r_uniform,edges,'FaceColor',color_uniform_scatter,'EdgeColor',color_uniform,'FaceAlpha',1,'Orientation','horizontal');
% xlim([0,20]);
% exportgraphics(gcf,'r_dist_RUM_uniform.eps','BackgroundColor','none','ContentType','vector');

%% fig. 5 - alppha uniform vs pareto (ESVT)
num_subjects = length(data(:,1));
figure;
scatter(alpha_uniform_ESVT,alpha_pareto_ESVT,160,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k','LineWidth',0.5);
hold on
x4 = 0:0.001:7;
y4=x4;
plot(x4,y4,'LineWidth',1.5,'LineStyle','--','Color','K'); 
hold on
errorbar(mean(alpha_uniform_ESVT),mean(alpha_pareto_ESVT),std(alpha_uniform_ESVT)./sqrt(num_subjects),'horizontal','CapSize',0,'LineWidth',3,'Color','k');
hold on
errorbar(mean(alpha_uniform_ESVT),mean(alpha_pareto_ESVT),std(alpha_pareto_ESVT)./sqrt(num_subjects),'vertical','CapSize',0,'LineWidth',3,'Color','k');
xlim([0,7]); ylim([0,7]);
% exportgraphics(gcf,'alpha_uniform_vs_pareto_esvt.eps','BackgroundColor','none','ContentType','vector');

% histogram inset
figure;
alpha_diff = alpha_uniform_ESVT - alpha_pareto_ESVT;
histogram(alpha_diff,20,'FaceColor',[0.7 0.7 0.7],'EdgeColor','k'); xlabel('\Delta \alpha'); ylabel('num. of subjects'); ylim([0 20]); xlim([-3 3]);
% exportgraphics(gcf,'alpha_diff_histogram.eps','BackgroundColor','none','ContentType','vector');
[p4,h4,stats4] = signrank(alpha_uniform_ESVT,alpha_pareto_ESVT,'tail','left');

figure;
scatter(subjects_r_uniform,subjects_r_pareto,160,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k','LineWidth',0.5);
hold on
x4 = 0:0.001:4;
y4=x4;
plot(x4,y4,'LineWidth',1.5,'LineStyle','--','Color','K'); 
hold on
errorbar(mean(subjects_r_uniform),mean(subjects_r_pareto),std(subjects_r_uniform)./sqrt(num_subjects),'horizontal','CapSize',0,'LineWidth',3,'Color','k');
hold on
errorbar(mean(subjects_r_uniform),mean(subjects_r_pareto),std(subjects_r_pareto)./sqrt(num_subjects),'vertical','CapSize',0,'LineWidth',3,'Color','k');
xlim([0,3.5]); ylim([0,3.5]);
% exportgraphics(gcf,'r_uniform_vs_pareto_rum.eps','BackgroundColor','none','ContentType','vector');

% histogram inset
figure;
r_diff = subjects_r_uniform - subjects_r_pareto;
histogram(r_diff,20,'FaceColor',[0.7 0.7 0.7],'EdgeColor','k'); xlabel('\Delta r'); ylabel('num. of subjects'); ylim([0 30]); xlim([-4 4]);
% exportgraphics(gcf,'r_diff_histogram.eps','BackgroundColor','none','ContentType','vector');
[p5,h5,stats5] = signrank(subjects_r_uniform,subjects_r_pareto,'tail','left');

%% fig. 6 - noise (violin plots)
noise_uniform_ESVT = data(:,5);
noise_uniform_RUM = data(:,14);
noise_pareto_ESVT = data(:,10);
noise_pareto_RUM = data(:,18);

figure;
noiseESVT_Table = array2table([noise_pareto_ESVT noise_uniform_ESVT ],...
                'VariableNames',{'Pareto','Uniform'});

v1=violinplot(noiseESVT_Table, [] ,'Width', 0.3,'ViolinColor',[color_pareto; color_uniform]);
xlim([0.5 2.5]); xticks([1 2]); set(gca,'FontSize',14); xticklabels({'Pareto','Uniform'});
ylabel('estimated noise (\epsilon)'); ylim([0 0.2]);
% exportgraphics(gcf,'noise_esvt_scatter.eps','BackgroundColor','none','ContentType','vector');
[p6,h6,stats6] = signrank(noise_uniform_ESVT,noise_pareto_ESVT);
%
figure;
noiseRUM_Table = array2table([noise_pareto_RUM noise_uniform_RUM],...
                'VariableNames',{'Pareto','Uniform'});

v1=violinplot(noiseRUM_Table, [] ,'Width', 0.3,'ViolinColor',[color_pareto; color_uniform]);
xlim([0.5 2.5]); xticks([1 2]); set(gca,'FontSize',14); xticklabels({'Pareto','Uniform'});
ylabel('estimated noise (\eta)'); ylim([0 1]);
% exportgraphics(gcf,'noise_rum.eps','BackgroundColor','none','ContentType','vector');
[p7,h7,stats7] = signrank(noise_uniform_RUM,noise_pareto_RUM);

%%
figure;
subplot(1,2,1);
scatter(noise_uniform_ESVT,noise_pareto_ESVT,160,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k','LineWidth',0.5);
hold on
x5 = 0:0.001:1;
y5=x5;
plot(x5,y5,'LineWidth',1.5,'LineStyle','--','Color','K'); 
xlim([0,0.2]); ylim([0,0.2]); set(gca,'FontSize',14);
xlabel('noise (\theta) - uniform (DN)','FontSize',14); ylabel('noise (\theta) - Pareto (DN)','FontSize',14);

subplot(1,2,2);
scatter(noise_uniform_RUM,noise_pareto_RUM,160,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k','LineWidth',0.5);
hold on
x5 = 0:0.001:2;
y5=x5;
plot(x5,y5,'LineWidth',1.5,'LineStyle','--','Color','K'); 
xlim([0,2]); ylim([0,2]);  set(gca,'FontSize',14);
xlabel('noise (\theta) - uniform (Power utility)','FontSize',14); ylabel('noise (\theta) - Pareto (Power utility)','FontSize',14);

