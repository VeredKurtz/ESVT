clear
clc 

rng('shuffle');

%%
color_pareto = [0 112 192]/255;
color_uniform = [191 32 37]/255;
color_pareto_scatter = [147 209 255]/255;
color_uniform_scatter = [238 156 158]/255;

%% Input dists
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
% figure;
% hist3(binary_unif_SV,[length(sv_array) length(sv_array)],'CdataMode','auto');
% xlabel('v1'); ylabel('v2'); 
% colorbar; view(2); clim([0 200]);
 % exportgraphics(gcf,'fig1_uniform_hist3.eps','BackgroundColor','none','ContentType','vector');
figure;
h = histogram2(binary_unif_SV(:,1),binary_unif_SV(:,2),sv_array,sv_array,'DisplayStyle','tile','ShowEmptyBins','on');
set(h, 'EdgeColor', 'none'); xlim([0 58]); ylim([0 58]); clim([0 200]); colormap(flipud(gray));
% exportgraphics(gcf,'fig1_uniform_hist3.eps','BackgroundColor','none','ContentType','vector');

% marginal dists
figure;
histogram(binary_unif_SV(:,1),40,'Normalization','probability');
xlabel('v1'); ylim([0 0.06]); 
 % exportgraphics(gcf,'fig1_uniform_marginal.eps','BackgroundColor','none','ContentType','vector');

% pareto dist
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
fig = clf();
ax1 = axes(fig); 
ax2 = copyobj(ax1,fig);

[N,xedges,yedges]=histcounts2(beta2(:,1),beta2(:,2));
hist_edges = 0:inc:max(x);
pc = histogram2(ax1,beta2(:,1),beta2(:,2),hist_edges,hist_edges,'DisplayStyle','tile','ShowEmptyBins','on');
set(pc, 'EdgeColor', 'none'); xlim([0 60]); ylim([0 60]); clim([0 200]); %colormap(flipud(gray));

fc = contour(ax2,yedges(1:end-1),xedges(1:end-1),N,5,'LineWidth',1); xlim([0 60]); ylim([0 60]); %colormap("abyss");
% Set colormaps
colormap(ax1,flipud(gray))
colormap(ax2,flipud(sky))
% Set all other properties of ax1 before moving on
% Finally, link the axis properties and turn off axis #2.
ax2.UserData = linkprop([ax1,ax2],...
     {'Position','InnerPosition','DataAspectRatio','xtick','ytick', ...
     'ydir','xdir','xlim','ylim'}); % add more props as needed
ax2.Visible = 'off';
% exportgraphics(gcf,'fig1_pareto_hist3.eps','BackgroundColor','none','ContentType','vector');


% marginal histogram
figure;
histogram(S(:,1),edges,'Normalization','probability');
xlabel('v1');
 % exportgraphics(gcf,'fig1_pareto_marginal.eps','BackgroundColor','none','ContentType','vector');

%% Neural encoding functions
% EU with r=1 (linear)
x = x';
v1_uniform_linear_encoding = x;
% EU with r=0.4
r1 = 0.4;
v1_uniform_concave_encoding = x(:,1).^r1;
% EU with r=2.5
r2 = 2.5;
v1_uniform_convex_encoding = x(:,1).^r2;
% ESVT with alpha=1, M = 15;
alpha1=1; M1= 15;
v1_uniform_esvt1 = (x(:,1).^alpha1)./(M1.^alpha1 + x(:,1).^alpha1);
% ESVT with alpha=4, M = 30;
alpha2=4; M2=30;
v1_uniform_esvt2 = (x(:,1).^alpha2)./(M2.^alpha2 + x(:,1).^alpha2);
% ESVT with alpha=7, M = 45;
alpha3=7; M3=45;
v1_uniform_esvt3 = (x(:,1).^alpha3)./(M3.^alpha3 + x(:,1).^alpha3); 



%%
figure;
plot1 = shadedErrorBar(x,v1_uniform_concave_encoding,0.05*max(x).^r1*ones(length(x),1), ...
                        'lineprops',{'-','Color','k','markerfacecolor',[0.7 0.7 0.7],'linewidth',2}); ylim([0 6]);
figure;
plot2 = shadedErrorBar(x,v1_uniform_linear_encoding,0.05*max(x)*ones(length(x),1), ...
                        'lineprops',{'-','Color','k','markerfacecolor',[0.7 0.7 0.7],'linewidth',2}); ylim([0 60]);
figure;
plot3 = shadedErrorBar(x,v1_uniform_convex_encoding,0.05*max(x)^r2*ones(length(x),1), ...
                        'lineprops',{'-','Color','k','markerfacecolor',[0.7 0.7 0.7],'linewidth',2}); ylim([0 30000]);
figure;
plot4 = shadedErrorBar(x,v1_uniform_esvt1,0.05*ones(length(x),1), ...
                        'lineprops',{'-','Color','k','markerfacecolor',[0.7 0.7 0.7],'linewidth',2}); ylim([0 1]);
figure;
plot5 = shadedErrorBar(x,v1_uniform_esvt2,0.05*ones(length(x),1), ...
                        'lineprops',{'-','Color','k','markerfacecolor',[0.7 0.7 0.7],'linewidth',2}); ylim([0 1]);
figure;
plot6 = shadedErrorBar(x,v1_uniform_esvt3,0.05*ones(length(x),1), ...
                        'lineprops',{'-','Color','k','markerfacecolor',[0.7 0.7 0.7],'linewidth',2}); ylim([0 1]);
% exportgraphics(gcf,'fig1_encoding_functions.eps','BackgroundColor','none','ContentType','vector');

%% choice functions

noise_EU_linear = normrnd(0,0.05*max(x),100000,4);
noise_EU_low_r = normrnd(0,0.05*(max(x).^r1),100000,4);
noise_EU_high_r = normrnd(0,0.05*(max(x).^r2),100000,4);
noise_ESVT_low = normrnd(0,0.05*(max(v1_uniform_esvt1)),100000,4);
noise_ESVT_med = normrnd(0,0.05*(max(v1_uniform_esvt2)),100000,4); 
noise_ESVT_high = normrnd(0,0.05*(max(v1_uniform_esvt3)),100000,4); 

num_trials = 1000;

for v1=1:40
    for v2=1:40
        current_v1 = sv_array(v1);
        current_v2 = sv_array(v2);
        
        % simulating 50-50 lotteries
        outcome1_lottery1 = unifrnd(0,2.*current_v1,num_trials,1);
        outcome2_lottery1 = 2.*current_v1-outcome1_lottery1;
        outcome1_lottery2 = unifrnd(0,current_v2,num_trials,1);
        outcome2_lottery2 = 2.*current_v2-outcome1_lottery2;

        for trial=1:num_trials

            % linear EU
            v1_linear = 0.5.*(outcome1_lottery1(trial,1) + noise_EU_linear(:,1)) + 0.5.*(outcome2_lottery1(trial,1) + noise_EU_linear(:,2));
            v2_linear = 0.5.*(outcome1_lottery2(trial,1) + noise_EU_linear(:,3)) + 0.5.*(outcome2_lottery2(trial,1) + noise_EU_linear(:,4));
            choice_func_eu_linear(v1,v2,trial) = length(find(v1_linear>=v2_linear))./length(v1_linear);

            % EU with r=0.4
            v1_low_r = 0.5.*(outcome1_lottery1(trial,1).^r1 + noise_EU_low_r(:,1)) + 0.5.*(outcome2_lottery1(trial,1).^r1 + noise_EU_low_r(:,2));
            v2_low_r = 0.5.*(outcome1_lottery2(trial,1).^r1 + noise_EU_low_r(:,3)) + 0.5.*(outcome2_lottery2(trial,1).^r1 + noise_EU_low_r(:,4));
            choice_func_eu_low_r(v1,v2,trial) = length(find(v1_low_r>=v2_low_r))./length(v1_low_r);

            % EU with r=2.5
            v1_high_r = 0.5.*(outcome1_lottery1(trial,1).^r2 + noise_EU_high_r(:,1)) + 0.5.*(outcome2_lottery1(trial,1).^r2 + noise_EU_high_r(:,2));
            v2_high_r = 0.5.*(outcome1_lottery2(trial,1).^r2 + noise_EU_high_r(:,3)) + 0.5.*(outcome2_lottery2(trial,1).^r2 + noise_EU_high_r(:,4));
            choice_func_eu_high_r(v1,v2,trial) = length(find(v1_high_r>=v2_high_r))./length(v1_high_r);

            % ESVT with M=15, alpha=1
            v1_esvt_low_alpha = 0.5*(((outcome1_lottery1(trial,1).^alpha1)./(M1.^alpha1 + outcome1_lottery1(trial,1).^alpha1)) + noise_ESVT_low(:,1)) + 0.5*(((outcome2_lottery1(trial,1).^alpha1)./(M1.^alpha1 + outcome2_lottery1(trial,1).^alpha1)) + noise_ESVT_low(:,2));
            v2_esvt_low_alpha = 0.5*(((outcome1_lottery2(trial,1).^alpha1)./(M1.^alpha1 + outcome1_lottery2(trial,1).^alpha1)) + noise_ESVT_low(:,3)) + 0.5*(((outcome2_lottery2(trial,1).^alpha1)./(M1.^alpha1 + outcome2_lottery2(trial,1).^alpha1)) + noise_ESVT_low(:,4));
            choice_func_esvt_low_alpha(v1,v2,trial) = length(find(v1_esvt_low_alpha>=v2_esvt_low_alpha))./length(v1_esvt_low_alpha);
    
            % ESVT with M=30, alpha=4
            v1_esvt_mid_alpha = 0.5*(((outcome1_lottery1(trial,1).^alpha2)./(M2.^alpha2 + outcome1_lottery1(trial,1).^alpha2)) + noise_ESVT_med(:,1)) + 0.5*(((outcome2_lottery1(trial,1).^alpha2)./(M2.^alpha2 + outcome2_lottery1(trial,1).^alpha2)) + noise_ESVT_med(:,2));
            v2_esvt_mid_alpha = 0.5*(((outcome1_lottery2(trial,1).^alpha2)./(M2.^alpha2 + outcome1_lottery2(trial,1).^alpha2)) + noise_ESVT_med(:,3)) + 0.5*(((outcome2_lottery2(trial,1).^alpha2)./(M2.^alpha2 + outcome2_lottery2(trial,1).^alpha2)) + noise_ESVT_med(:,4));
            choice_func_esvt_mid_alpha(v1,v2,trial) = length(find(v1_esvt_mid_alpha>=v2_esvt_mid_alpha))./length(v1_esvt_mid_alpha);

            % ESVT with M=45, alpha=7
            v1_esvt_high_alpha = 0.5*(((outcome1_lottery1(trial,1).^alpha3)./(M3.^alpha3 + outcome1_lottery1(trial,1).^alpha3)) + noise_ESVT_high(:,1)) + 0.5*(((outcome2_lottery1(trial,1).^alpha3)./(M3.^alpha3 + outcome2_lottery1(trial,1).^alpha3)) + noise_ESVT_high(:,2));
            v2_esvt_high_alpha = 0.5*(((outcome1_lottery2(trial,1).^alpha3)./(M3.^alpha3 + outcome1_lottery2(trial,1).^alpha3)) + noise_ESVT_high(:,3)) + 0.5*(((outcome2_lottery2(trial,1).^alpha3)./(M3.^alpha3 + outcome2_lottery2(trial,1).^alpha3)) + noise_ESVT_high(:,4));
            choice_func_esvt_high_alpha(v1,v2,trial) = length(find(v1_esvt_high_alpha>=v2_esvt_high_alpha))./length(v1_esvt_high_alpha);

            clear v1_esvt_high_alpha v2_esvt_high_alpha v1_esvt_mid_alpha v2_esvt_mid_alpha v1_esvt_low_alpha v2_esvt_low_alpha
            clear v1_linear v2_linear v1_low_r v2_low_r v1_high_r v2_high_r
        
        end
        clear outcome1_lottery1 outcome2_lottery1 outcome1_lottery2 outcome2_lottery2

    end
end

%% choice function figure

% means
mean_choice_func_eu_low_r = squeeze(nanmean(choice_func_eu_low_r,3));
mean_choice_func_eu_linear = squeeze(nanmean(choice_func_eu_linear,3));
mean_choice_func_eu_high_r = squeeze(nanmean(choice_func_eu_high_r,3));
mean_choice_func_esvt_low_alpha = squeeze(nanmean(choice_func_esvt_low_alpha,3));
mean_choice_func_esvt_mid_alpha = squeeze(nanmean(choice_func_esvt_mid_alpha,3));
mean_choice_func_esvt_high_alpha = squeeze(nanmean(choice_func_esvt_high_alpha,3));

%
figure;
subplot(3,2,1);
imagesc(sv_array,sv_array,mean_choice_func_eu_low_r);
set(gca,'YDir','normal'); set(gca,'XAxisLocation','bottom');
colormap("gray"); 
subplot(3,2,3);
imagesc(sv_array,sv_array,mean_choice_func_eu_linear);
set(gca,'YDir','normal'); set(gca,'XAxisLocation','bottom');
colormap("gray"); 
subplot(3,2,5);
imagesc(sv_array,sv_array,mean_choice_func_eu_high_r);
set(gca,'YDir','normal'); set(gca,'XAxisLocation','bottom');
colormap("gray"); 

subplot(3,2,2);
imagesc(sv_array,sv_array,mean_choice_func_esvt_low_alpha);
set(gca,'YDir','normal'); set(gca,'XAxisLocation','bottom');
colormap("gray"); 
subplot(3,2,4);
imagesc(sv_array,sv_array,mean_choice_func_esvt_mid_alpha);
set(gca,'YDir','normal'); set(gca,'XAxisLocation','bottom');
colormap("gray"); 
subplot(3,2,6);
imagesc(sv_array,sv_array,mean_choice_func_esvt_high_alpha);
set(gca,'YDir','normal'); set(gca,'XAxisLocation','bottom'); 
colormap("gray"); 
% exportgraphics(gcf,'fig1_choice_functions.eps','BackgroundColor','none','ContentType','vector');
