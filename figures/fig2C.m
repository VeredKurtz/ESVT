clear
clc 

rng('shuffle');

%%
color_pareto = [0 112 192]/255;
color_uniform = [191 32 37]/255;
color_pareto_scatter = [147 209 255]/255;
color_uniform_scatter = [238 156 158]/255;

%% RN
x = 0:0.001:60;
num_sims = 100000;
choice_set_size = 2;

% uniform dist
sv_inc = max(x)./39;
sv_array = 0:sv_inc:max(x);
sv_array = sv_array';
repetitions = num_sims./length(sv_array);
rep_sv = repelem(sv_array,repetitions);
for i=1:choice_set_size
    random_order_sv(:,i) = randperm(length(rep_sv));
    binary_unif_SV(:,i) = rep_sv(random_order_sv(:,i));
end

% bivariate hists
figure;
h = histogram2(binary_unif_SV(:,1),binary_unif_SV(:,2),sv_array,sv_array,'DisplayStyle','tile','ShowEmptyBins','on');
set(h, 'EdgeColor', 'none'); xlim([0 58]); ylim([0 58]); clim([0 200]); colormap(flipud(gray));
% exportgraphics(gcf,'fig1_uniform_hist3.eps','BackgroundColor','none','ContentType','vector');

% marginal dists
figure;
histogram(binary_unif_SV(:,1),40,'Normalization','probability');
xlabel('v1'); ylim([0 0.06]); 
 % exportgraphics(gcf,'fig1_uniform_marginal.eps','BackgroundColor','none','ContentType','vector');

%% pareto dist
beta = 3; k=1/beta; theta = 0; 
a = 1; % for the gamma distribution
b = 1; % for the gamma distribution
sigma = mean(x)./(gamma((beta-1)./beta).*gamma((beta+1)./beta)); 

for trial=1:num_sims
    Z = gamrnd(a,b);
    for i=1:choice_set_size
        U(i) = exprnd(1);
        S(trial,i)= theta + sigma * (U(i)/Z)^(1/beta);
    end
end

inc = max(x)./39;
edges = 0:inc:max(x);
    
% bivariate histogram
figure;
beta2 = [S(:,1) S(:,2)];
inc = max(x)./39;
hist3(beta2,'edges',{0:inc:max(x) 0:inc:max(x)},'CdataMode','auto');
view(2); clim([0 200]); xlim([0 60]); ylim([0 60]); 
%colorbar;
% exportgraphics(gcf,'fig1_pareto_hist3.eps','BackgroundColor','none','ContentType','vector');

% contour plot
figure;
[N,xedges,yedges]=histcounts2(beta2(:,1),beta2(:,2));
hist_edges = 0:inc:max(x);
pc = histogram2(beta2(:,1),beta2(:,2),hist_edges,hist_edges,'DisplayStyle','tile','ShowEmptyBins','on');
set(pc, 'EdgeColor', 'none'); xlim([0 60]); ylim([0 60]); clim([0 250]); %colormap(flipud(gray));
colormap(flipud(gray));
% exportgraphics(gcf,'fig3_pareto_hist3_RN.eps','BackgroundColor','none','ContentType','vector')

% marginal histogram
figure;
histogram(S(:,1),edges,'Normalization','probability');
xlabel('v1'); ylim([0 0.07]);
exportgraphics(gcf,'fig3_pareto_marginal_RN.eps','BackgroundColor','none','ContentType','vector');

%% RS
rho = 1.61;
x2 = x.^rho;
num_sims = 100000;
choice_set_size = 2;

% uniform dist
sv_inc = max(x2)./39;
sv_array = 0:sv_inc:max(x2);
sv_array = sv_array';
repetitions = num_sims./length(sv_array);
rep_sv = repelem(sv_array,repetitions);
for i=1:choice_set_size
    random_order_sv(:,i) = randperm(length(rep_sv));
    binary_unif_SV(:,i) = rep_sv(random_order_sv(:,i));
end

% bivariate hists
figure;
h = histogram2(binary_unif_SV(:,1),binary_unif_SV(:,2),sv_array,sv_array,'DisplayStyle','tile','ShowEmptyBins','on');
set(h, 'EdgeColor', 'none'); xlim([0 58]); ylim([0 58]); clim([0 200]); colormap(flipud(gray));
% exportgraphics(gcf,'fig1_uniform_hist3.eps','BackgroundColor','none','ContentType','vector');

% marginal dists
figure;
histogram(binary_unif_SV(:,1),40,'Normalization','probability');
xlabel('v1'); ylim([0 0.06]); 
% exportgraphics(gcf,'fig1_uniform_marginal.eps','BackgroundColor','none','ContentType','vector');

%% pareto dist
beta = 3; k=1/beta; theta = 0; 
a = 1; % for the gamma distribution
b = 1; % for the gamma distribution
sigma = mean(x2)./(gamma((beta-1)./beta).*gamma((beta+1)./beta)); 

for trial=1:num_sims
    Z = gamrnd(a,b);
    for i=1:choice_set_size
        U(i) = exprnd(1);
        S(trial,i)= theta + sigma * (U(i)/Z)^(1/beta);
    end
end

inc = max(x2)./39;
edges = 0:inc:max(x2);
    
% bivariate histogram
figure;
beta2 = [S(:,1) S(:,2)];
inc = max(x2)./39;
hist3(beta2,'edges',{0:inc:max(x2) 0:inc:max(x2)},'CdataMode','auto');
view(2); clim([0 200]); xlim([0 max(x2)]); ylim([0 max(x2)]); 
%colorbar;
% exportgraphics(gcf,'fig3_pareto_hist3.eps','BackgroundColor','none','ContentType','vector');

% contour plot
figure;
[N,xedges,yedges]=histcounts2(beta2(:,1),beta2(:,2));
hist_edges = 0:inc:max(x2);
pc = histogram2(beta2(:,1),beta2(:,2),hist_edges,hist_edges,'DisplayStyle','tile','ShowEmptyBins','on');
set(pc, 'EdgeColor', 'none'); xlim([0 max(x2)]); ylim([0 max(x2)]); clim([0 250]);
colormap(flipud(gray));
% exportgraphics(gcf,'fig3_pareto_hist3_RS.eps','BackgroundColor','none','ContentType','vector')

% marginal histogram
figure;
histogram(S(:,1),edges,'Normalization','probability');
xlabel('v1'); ylim([0 0.07])
exportgraphics(gcf,'fig3_pareto_marginal_RS.eps','BackgroundColor','none','ContentType','vector');

%% RA
rho = 0.43;
x2 = x.^rho;
num_sims = 100000;
choice_set_size = 2;

% uniform dist
sv_inc = max(x2)./39;
sv_array = 0:sv_inc:max(x2);
sv_array = sv_array';
repetitions = num_sims./length(sv_array);
rep_sv = repelem(sv_array,repetitions);
for i=1:choice_set_size
    random_order_sv(:,i) = randperm(length(rep_sv));
    binary_unif_SV(:,i) = rep_sv(random_order_sv(:,i));
end

% bivariate hists
figure;
h = histogram2(binary_unif_SV(:,1),binary_unif_SV(:,2),sv_array,sv_array,'DisplayStyle','tile','ShowEmptyBins','on');
set(h, 'EdgeColor', 'none'); xlim([0 58]); ylim([0 58]); clim([0 200]); colormap(flipud(gray));
% exportgraphics(gcf,'fig1_uniform_hist3.eps','BackgroundColor','none','ContentType','vector');

% marginal dists
figure;
histogram(binary_unif_SV(:,1),40,'Normalization','probability');
xlabel('v1'); ylim([0 0.06]); 
 % exportgraphics(gcf,'fig1_uniform_marginal_RS.eps','BackgroundColor','none','ContentType','vector');

%% pareto dist
beta = 3; k=1/beta; theta = 0; 
a = 1; % for the gamma distribution
b = 1; % for the gamma distribution
sigma = mean(x2)./(gamma((beta-1)./beta).*gamma((beta+1)./beta)); 

for trial=1:num_sims
    Z = gamrnd(a,b);
    for i=1:choice_set_size
        U(i) = exprnd(1);
        S(trial,i)= theta + sigma * (U(i)/Z)^(1/beta);
    end
end

inc = max(x2)./39;
edges = 0:inc:max(x2);
    
% bivariate histogram
figure;
beta2 = [S(:,1) S(:,2)];
inc = max(x2)./39;
hist3(beta2,'edges',{0:inc:max(x2) 0:inc:max(x2)},'CdataMode','auto');
view(2); clim([0 200]); xlim([0 max(x2)]); ylim([0 max(x2)]); 
%colorbar;
% exportgraphics(gcf,'fig1_pareto_hist3.eps','BackgroundColor','none','ContentType','vector');

% contour plot
figure;
[N,xedges,yedges]=histcounts2(beta2(:,1),beta2(:,2));
hist_edges = 0:inc:max(x2);
pc = histogram2(beta2(:,1),beta2(:,2),hist_edges,hist_edges,'DisplayStyle','tile','ShowEmptyBins','on');
set(pc, 'EdgeColor', 'none'); xlim([0 max(x2)]); ylim([0 max(x2)]); clim([0 250]);
colormap(flipud(gray));
% exportgraphics(gcf,'fig3_pareto_hist3_RA.eps','BackgroundColor','none','ContentType','vector')

% marginal histogram
figure;
histogram(S(:,1),edges,'Normalization','probability');
xlabel('v1'); ylim([0 0.07]);
exportgraphics(gcf,'fig3_pareto_marginal_RA.eps','BackgroundColor','none','ContentType','vector');
