function [sv_binary_exp,sv_exp_sixOpt] = pareto_dist_generator_subject_specific_newVersion(subjects_alpha,means_of_distributions,max_SV_value)

% DRAWS FOR THE PARETO DISTRIBUTIONS IN THE ESVT EXPERIMENT
%  In this version cast the joint Pareto dist. based on Gamma and Exp distributions.
%  I use each subject's elicited alpha parameter,
%  and based on the SV distributions in the Uniform case,
%  I build each subject's Pareto distribution.
%  I first create a million draws distribution, calculate its std,
%  and then simulate this distribution in a smaller 380 drwas distribution
%  for the experiment.

num_of_subjects = length(subjects_alpha);

%  According to Bucher & Brandenburger (2021) 
%  the equation of E(S) at p. 24 in the Appendix Derivation of Moments,
%  once we fix E(S) to equal each subject's mean from the uniform dist., 
%  like in our experiment, we can fix mu=0, and solve for sigma. 

%  beta=3
gamma1 = gamma(2/3);
gamma2 = gamma(1/3);
gamma3 = gamma(5/3);
correl = (gamma1.^(-2) - gamma2.^(-1))./gamma3;
sig_beta3 = 17.5./(gamma(2/3).*gamma(4/3));
%  beta=2
gamma11 = gamma(1/2);
gamma21 = gamma(0);
gamma31 = gamma(2);
correl11 = (gamma11.^(-2) - gamma21.^(-1))./gamma31;

%% parameters
beta = 3; % decided on ESVT group meeting on 25 May 2021
k=1/beta;
theta = 0; % we set mu (theta) to be 0
a = 1; % for the gamma distribution
b = 1; % for the gamma distribution

% calculate each subject's sigma, following Bucher and Brandenburger (2021)
for sub=1:num_of_subjects
    sigma(sub) = means_of_distributions(sub)./(gamma((beta-1)./beta).*gamma((beta+1)./beta)); 
end
    
%% BINARY SETS

for sub=1:num_of_subjects
    for trial=1:100000
        Z = gamrnd(a,b);
        for i=1:2
            U(i) = exprnd(1);
            S(trial,i,sub)= theta + sigma(sub) * (U(i)/Z)^(1/beta);
        end
    end
    
    % calculate how much we trancuate from the distribution
    num_of_obs_left_out_s2 = sum(S(:,2,sub)>max_SV_value(sub));
    trancuatedProb_s2(sub) = num_of_obs_left_out_s2./length(S(:,2,sub));
        
    num_of_obs_left_out_s1 = sum(S(:,1,sub)>max_SV_value(sub));
    trancuatedProb_s1(sub) = num_of_obs_left_out_s1./length(S(:,1,sub));
    
    % descriptive statistics of distributions
    mean_s2_binary(sub) = mean(S(:,2,sub));
    std_s2_binary(sub) = std(S(:,2,sub));
    
    mean_s1_binary(sub) = mean(S(:,1,sub));
    std_s1_binary(sub) = std(S(:,1,sub));
    
    pearson(sub) = corr(S(:,1,sub),S(:,2,sub));

    inc = max_SV_value(sub)./100;
    edges = 0:inc:max_SV_value(sub);
    
%     subplot(2,7,sub);
%     histogram(S(2,:,sub),edges,'Normalization','probability');
%     title({['dist of S2, sub ' num2str(sub) ',beta=' num2str(beta)], ...
%         ['mu=0,sigma=' num2str(sigma(sub))], ...
%         ['mean=' num2str(mean_s2_binary(sub))], ...
%         ['std=' num2str(std_s2_binary(sub))],...
%         ['median=' num2str(median(S(:,2,sub)))], ...
%         ['trancuated prob =' num2str(trancuatedProb_s2(sub))]},'FontSize',7);
%         ylim([0 0.025]);
    
%     subplot(2,7,sub+7);
%     histogram(S(1,:,sub),edges,'Normalization','probability');
%     title({['dist of S1, sub ' num2str(sub) ',beta=' num2str(beta)],...
%         ['mu=0,sigma=' num2str(sigma(sub))], ...
%         ['mean=' num2str(mean_s1_binary(sub))], ...
%         ['std=' num2str(std_s1_binary(sub))],...
%         ['median='  num2str(median(S(:,1,sub)))], ...
%         ['trancuated prob =' num2str(trancuatedProb_s1(sub))]},'FontSize',7);
%         ylim([0 0.025]);
    
end

% sgtitle('marginal distributions');

%% Bivariate bistogram
% figure()
% for sub = 1:num_of_subjects
%     beta2 = [S(1,:,sub)' S(2,:,sub)'];
%     inc = max_SV_value(sub)./50;
%     subplot(2,4,sub);
%     hist3(beta2,'edges',{0:inc:max_SV_value(sub) 0:inc:max_SV_value(sub)},'CdataMode','auto');
%     view(2)
%     colorbar;
%     title(['sub ' num2str(sub) ',corr= ' num2str(round(pearson(sub),3))]); xlabel('s1'); ylabel('s2');
% end
% sgtitle('bivariate distributions');

%% Create a small 600-pairs of SVs for the experiment
% The small set should deviate by 2% from the mean and std of the large
% distribution

deviation_threshold = 0.02; % in percentile
% We will need only 320 lotteries in the experiment. 
% We cast extra values because we willneed to trancuate some SVs
num_of_draws = 600; 
sv_binary_exp = zeros(num_of_draws,2,num_of_subjects);
% start a subjects loop
for sub=1:num_of_subjects
    
    threshold_mean_s2 = deviation_threshold.*mean_s2_binary(sub);
    threshold_std_s2 = deviation_threshold.*mean_s2_binary(sub);
    
    threshold_mean_s1 = deviation_threshold.*mean_s1_binary(sub);
    threshold_std_s1 = deviation_threshold.*mean_s1_binary(sub);
    
    deviation_mean_s2 = 10000; % set arbitrary values
    deviation_std_s2 = 10000; % set arbitrary values
    deviation_mean_s1 = 10000; % set arbitrary values
    deviation_std_s1 = 10000; % set arbitrary values
    
    % draw the small 600 lotteries dist 
    while deviation_mean_s2>threshold_mean_s2 || deviation_std_s2>threshold_std_s2 || deviation_mean_s1>threshold_mean_s1 || deviation_std_s1>threshold_std_s1
        Z = gamrnd(a,b,num_of_draws,1);
        U = exprnd(1,num_of_draws,2);
        for i=1:2
            sv_binary_exp(:,i,sub)= theta + sigma(sub) .* (U(:,i)./Z(:,1)).^(1/beta);
        end
        meanS2(sub) = mean(sv_binary_exp(:,2,sub));
        stdS2(sub) = std(sv_binary_exp(:,2,sub));
        deviation_mean_s2 = abs(meanS2(sub)-mean_s2_binary(sub)); % compute deviation 
        deviation_std_s2 = abs(stdS2(sub)-std_s2_binary(sub)); % compute deviation 
        meanS1(sub) = mean(sv_binary_exp(:,1,sub));
        stdS1(sub) = std(sv_binary_exp(:,1,sub));
        deviation_mean_s1 = abs(meanS1(sub)-mean_s1_binary(sub)); % compute deviation 
        deviation_std_s1 = abs(stdS1(sub)-std_s1_binary(sub)); % compute deviation 
    end
    
    % create a table for export
    s1_exp_table = array2table(sv_binary_exp(:,1,sub),'VariableNames',{'s1'});
    s2_exp_table = array2table(sv_binary_exp(:,2,sub),'VariableNames',{'s2'});

end


%% SIX OPTIONS SETS 
% Here we create the distributions for the six-options Pareto lotteries
% The procedure is almost identical to the binary case,
% only we expand the draws from the quantile function to six options

for sub = 1:num_of_subjects
    Z = gamrnd(a,b,100000,1);
    U = exprnd(1,100000,6);
    for i=1:6
        sixOptionsSV(:,i,sub) = theta + sigma(sub) .* (U(:,i)./Z(:,1)).^(1/beta);
    end 
end


%% Calculate means and stds of the sextuplet distributions
mean_sextuplet = mean(sixOptionsSV,1);
std_sextuplet = std(sixOptionsSV,0,1);
for sub=1:num_of_subjects
    correlations_sextuplet(:,:,sub) = corr(sixOptionsSV(:,:,sub));
end

%% Create five small 600-sextuplets of SVs for the experiment
% The small sets should deviate by 2% from the mean and std of the large
% distribution

deviation_threshold = 0.02;
num_of_draws = 600;
sv_exp_sixOpt = zeros(num_of_draws,6,num_of_subjects);

deviation_mean = 10000.*ones(num_of_subjects,6); % set arbitrary values
deviation_std = 10000.*ones(num_of_subjects,6); % set arbitrary values

% start a subjects loop
for sub=1:num_of_subjects
    threshold_mean_sixOpt = deviation_threshold.*mean_sextuplet(:,:,sub);
    threshold_std_sixOpt = deviation_threshold.*std_sextuplet(:,:,sub);

    Z = gamrnd(a,b,num_of_draws,1);
    for j=1:6
    % next draw the small 600 lotteries dist for s1-s6
        while deviation_mean(sub,j)>threshold_mean_sixOpt(1,j) || deviation_std(sub,1)>threshold_std_sixOpt(1,j) 
            U = exprnd(1,num_of_draws,1);
            sv_exp_sixOpt(:,j,sub)= theta + sigma(sub) .* (U(:,1)./Z(:,1)).^(1/beta);
            means_exp_sixOpt(sub,j) = mean(sv_exp_sixOpt(:,j,sub));
            stds_exp_sixOpt(sub,j) = std(sv_exp_sixOpt(:,j,sub));

            deviation_mean(sub,j) = abs((means_exp_sixOpt(sub,j)-mean_sextuplet(1,j,sub))); % compute deviation
            deviation_std(sub,j) = abs((stds_exp_sixOpt(sub,j)-std_sextuplet(1,j,sub))); % compute deviation
        end     
    end
end

