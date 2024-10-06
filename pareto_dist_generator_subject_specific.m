function [sv_binary_exp,sv_exp_sixOpt] = pareto_dist_generator_subject_specific(subjects_alpha,means_of_distributions,max_SV_value)

% 
% subjects_alpha = [0.1 0.3 1.1];
% means_of_distributions = [1.25298292730746;2.20771494618988;45.6789756384477];
% max_SV_value = [1.50596585461492 3.41542989237976 90.3579512768954];

% DRAWS FOR THE PARETO DISTRIBUTIONS IN THE ESVT EXPERIMENT
%  In this version I don't use the Matlab built-in function for Pareto
%  dist, and draw from a distribution I created, for both s1 and s2 in the
%  binary case, and for si, i=1,...,6 in the six-options set.
%  I use each subject's elicited alpha parameter,
%  and based on the SV distributions in the Uniform case,
% I build each subject's Pareto distribution.
% I first create a million draws distribution, calculate its std,
% and then simulate this distribution in a smaller 320 drwas distribution
% for the experiment.

num_of_subjects = length(subjects_alpha);

% According to Bucher & Brandenburger (2021) 
% the equation of E(S) at p. 24 in the Appendix Derivation of Moments,
% once we fix E(S) to equal each subject's mean from the uniform dist., 
% like in our experiment, we can fix mu=0, and solve for sigma. 
% beta=3
gamma1 = gamma(2/3);
gamma2 = gamma(1/3);
gamma3 = gamma(5/3);
corr = (gamma1.^(-2) - gamma2.^(-1))./gamma3;
sig_beta3 = 17.5./(gamma(2/3).*gamma(4/3));
% beta=2
gamma11 = gamma(1/2);
gamma21 = gamma(0);
gamma31 = gamma(2);
corr11 = (gamma11.^(-2) - gamma21.^(-1))./gamma31;

%% parameters
beta = 3; % decided on ESVT group meeting on 25 May 2021
k=1/beta;
theta = 0; % we set mu (theta) to be 0

% calculate each subject's sigma, following Bucher and Brandenburger (2021)
for sub=1:num_of_subjects
    sigma(sub) = means_of_distributions(sub)./(gamma((beta-1)./beta).*gamma((beta+1)./beta)); 
end
    

%% BINARY SETS: solve for s2
% according to Bucher & Brandenburger (2021), s2 is taken from a Pareto
% type III dist (equation 7). We can isolate s2, and get - 
% s2 = (1/F_s2 - 1).^(-1/beta).*sigma_j+mu_j
% we set mu = 0 and sigma is subject-specific, so - 
% s2 = sigma.*(1/F_s2 - 1).^(-1/beta)
% where F_s2 are randomly drawn probabilities

for sub = 1:num_of_subjects
    % implementation
    % draw random values for the marginal cdf
    % I cast 1,000,000 samples to capture how the distribution 
    % looks like, so that when I cast smaller 300-lotteries sets, 
    % i could capture the moments of the distribution
    cdf_marginal(:,sub) = unifrnd(0,1,[100000,1]);

    % now solve for s2
    for i=1:length(cdf_marginal)
        s2(i,sub) = sigma(sub).*((1./cdf_marginal(i,sub))-1).^(1./beta);
    end
    
    % sample s1
    % Based on Bucher & Brandenburger (2021), given the conditional CDF,
    % and given that mu_j=mu_i=0
    % s1 is defined by - 
    % s1 = sigma_i[((1-Fs2|s1).^(-0.5) - 1].^(1/beta)[1+(s2./sigma_j).^beta].^(1/beta) 
    % The conditional cdf is also a random draw of values between 0 and 1.
    % all the other expressions in s1 are parameter or the draws from s2

    % draw random values for the conditional cdf
    cdf_conditional(:,sub) = unifrnd(0,1,[100000,1]);

    % now solve for s1
    for i=1:length(s2)
        s1(i,sub) = sigma(sub).*((((1-cdf_conditional(i,sub)).^(-0.5)-1).^k).*(1+(s2(i,sub)./sigma(sub)).^beta).^k);
    end

    % calculate how much we trancuate from the distribution
    num_of_obs_left_out_s2 = sum(s2(:,sub)>max_SV_value(sub));
    trancuatedProb_s2(sub) = num_of_obs_left_out_s2./length(s2(:,sub));
        
    num_of_obs_left_out_s1 = sum(s1(:,sub)>max_SV_value(sub));
    trancuatedProb_s1(sub) = num_of_obs_left_out_s1./length(s1(:,sub));
    
    % descriptive statistics of distributions
    mean_s2_binary(sub) = mean(s2(:,sub));
    std_s2_binary(sub) = std(s2(:,sub));

    inc = max_SV_value(sub)./100;
    edges = 0:inc:max_SV_value(sub);
    
%     subplot(2,7,sub);
%     histogram(s2(:,sub),edges,'Normalization','probability');
%     title({['dist of S2, sub ' num2str(sub) ',beta=' num2str(beta)], ...
%         ['mu=0,sigma=' num2str(sigma(sub))], ...
%         ['mean=' num2str(mean_s2_binary(sub))], ...
%         ['std=' num2str(std_s2_binary(sub))],...
%         ['median=' num2str(median(s2(:,sub)))], ...
%         ['trancuated prob =' num2str(trancuatedProb_s2(sub))]},'FontSize',7);
% %        ylim([0 0.04]);
    
    mean_s1_binary(sub) = mean(s1(:,sub));
    std_s1_binary(sub) = std(s1(:,sub));

%     subplot(2,7,sub+7);
%     histogram(s1(:,sub),edges,'Normalization','probability');
%     title({['dist of S1, sub ' num2str(sub) ',beta=' num2str(beta)],...
%         ['mu=0,sigma=' num2str(sigma(sub))], ...
%         ['mean=' num2str(mean_s1_binary(sub))], ...
%         ['std=' num2str(std_s1_binary(sub))],...
%         ['median='  num2str(median(s2(:,sub)))], ...
%         ['trancuated prob =' num2str(trancuatedProb_s1(sub))]},'FontSize',7);
% %      ylim([0 0.04]);
    

end
% sgtitle('marginal distributions');

%% Bivariate bistogram
% figure()
% for sub = 1:num_of_subjects
%     beta2 = [s1(:,sub) s2(:,sub)];
%     inc = max_SV_value(sub)./50;
%     subplot(2,4,sub);
%     hist3(beta2,'edges',{0:inc:max_SV_value(sub) 0:inc:max_SV_value(sub)},'CdataMode','auto');
%     view(2)
%     colorbar;
%     title(['sub ' num2str(sub)]); xlabel('s1'); ylabel('s2');
% end
% sgtitle('bivariate distributions');

%% Create a small 380-pairs of SVs for the experiment
% The small set should deviate by 1% from the mean and std of the large
% distribution

deviation_threshold = 0.01; % in percentile
% We will need only 320 lotteries in the experiment. 
% We cast extra values because we willneed to trancuate some SVs
num_of_draws = 380; 
sv_binary_exp = zeros(num_of_draws,2,num_of_subjects);
% start a subjects loop
for sub=1:num_of_subjects
    
    threshold_mean_s2 = deviation_threshold.*mean_s2_binary(sub);
    threshold_std_s2 = deviation_threshold.*mean_s2_binary(sub);
    
    threshold_mean_s1 = deviation_threshold.*mean_s1_binary(sub);
    threshold_std_s1 = deviation_threshold.*mean_s1_binary(sub);
    
    deviation_mean_s2 = 10000; % set arbitrary values
    deviation_std_s2 = 10000; % set arbitrary values
    % first draw the small 320 lotteries dist for s2
    while deviation_mean_s2>threshold_mean_s2 || deviation_std_s2>threshold_std_s2
        cdf_marginal_exp(:,sub) = unifrnd(0,1,[num_of_draws,1]);
        for i=1:num_of_draws
            sv_binary_exp(i,2,sub) = sigma(sub).*((1./cdf_marginal_exp(i,sub))-1).^(1./beta);
        end
        meanS2(sub) = mean(sv_binary_exp(:,2,sub));
        stdS2(sub) = std(sv_binary_exp(:,2,sub));
        deviation_mean_s2 = abs(meanS2(sub)-mean_s2_binary(sub)); % compute deviation 
        deviation_std_s2 = abs(stdS2(sub)-std_s2_binary(sub)); % compute deviation 
    end
    
    % next draw the small 320 lotteries dist for s1
    deviation_mean_s1 = 10000; % set arbitrary values
    deviation_std_s1 = 10000; % set arbitrary values
    % first draw the small 320 lotteries dist for s1
    while deviation_mean_s1>threshold_mean_s1 || deviation_std_s1>threshold_std_s1
        cdf_conditional_exp(:,sub) = unifrnd(0,1,[num_of_draws,1]);
        for i=1:num_of_draws
            sv_binary_exp(i,1,sub) = sigma(sub).*((((1-cdf_conditional_exp(i,sub)).^(-0.5)-1).^k).*(1+(sv_binary_exp(i,2,sub)./sigma(sub)).^beta).^k);
        end
        meanS1(sub) = mean(sv_binary_exp(:,1,sub));
        stdS1(sub) = std(sv_binary_exp(:,1,sub));
        deviation_mean_s1 = abs(meanS1(sub)-mean_s1_binary(sub)); % compute deviation 
        deviation_std_s1 = abs(stdS1(sub)-std_s1_binary(sub)); % compute deviation 
    end   
      
    % create a table for export
    s1_exp_table = array2table(sv_binary_exp(:,1,sub),'VariableNames',{'s1'});
    s2_exp_table = array2table(sv_binary_exp(:,2,sub),'VariableNames',{'s2'});
    output = [s1_exp_table, s2_exp_table];
%     fileName = ['sv_distributions_Pareto_binary_sets_sub' num2str(sub) '.csv'];
%     writetable(output,fileName);
    
    clear output

end


%% SIX OPTIONS SETS 
% Here we create the distributions for the six-options Pareto lotteries
% The procedure is almost identical to the binary case,
% only we expand the draws from the quantile function to six options

for sub = 1:num_of_subjects
    % implementation
    % draw random values for the marginal cdf with a million draws
    cdf_sixOptionsSV(:,1,sub) = unifrnd(0,1,[1000000,1]);

    % now solve for s1
    for i=1:length(cdf_sixOptionsSV)
        sixOptionsSV(i,1,sub) = sigma(sub).*((1./cdf_sixOptionsSV(i,1,sub))-1).^(1./beta);
    end
%     subplot(7,6,1+(sub-1).*6);
    inc = max_SV_value(sub)./50;
    edges = 0:inc:max_SV_value(sub);
%     histogram(sixOptionsSV(:,1,sub),edges,'Normalization','probability');
%     title(['sub' num2str(sub)],'FontSize',10);

    % now start a loop for the other five draws
    % where sj = sigma_j[(1-CDF_j)^(-1/j) -
    % 1][(1+sum(si-mu_i/sigma_i)^beta)]^(1/beta) + mu_j  (note tha mu+i=mu_j=0)
    for j=2:6
       cdf_sixOptionsSV(:,j,sub) = unifrnd(0,1,[1000000,1]);
       for i=1:length(cdf_sixOptionsSV)
           summationOfSVs=0;
           for z=1:j-1
               summationOfSVs = summationOfSVs+(sixOptionsSV(i,z,sub)./sigma(sub)).^beta;
           end
           sixOptionsSV(i,j,sub) = sigma(sub).*((((1-cdf_sixOptionsSV(i,j,sub)).^(-1./j)-1).^k).*(1+summationOfSVs).^k);
       end
%        subplot(7,6,j+(sub-1).*6);
%        inc = max_SV_value(sub)./50;
%        edges = 0:inc:max_SV_value(sub);
%        histogram(sixOptionsSV(:,j,sub),edges,'Normalization','probability');
    end
end


%% Calculate means and stds of the sextuplet distributions

mean_sextuplet = mean(sixOptionsSV,1);
std_sextuplet = std(sixOptionsSV,0,1);
for sub=1:num_of_subjects
    correlations_sextuplet(:,:,sub) = corr(sixOptionsSV(:,:,sub));
end

%% Create five small 380-sextuplets of SVs for the experiment
% The small sets should deviate by 1% from the mean and std of the large
% distribution

deviation_threshold = 0.01;
num_of_draws = 450;
sv_exp_sixOpt = zeros(num_of_draws,6,num_of_subjects);

% start a subjects loop
for sub=1:num_of_subjects
      
    threshold_mean_sixOpt = deviation_threshold.*mean_sextuplet(:,:,sub);
    threshold_std_sixOpt = deviation_threshold.*std_sextuplet(:,:,sub);
    
    deviation_mean_s1 = 10000; % set arbitrary values
    deviation_std_s1 = 10000; % set arbitrary values
    % first draw the small 320 lotteries dist for s1
    while deviation_mean_s1>threshold_mean_sixOpt(1,1) || deviation_std_s1>threshold_std_sixOpt(1,1)
        cdf_exp_sixOpt(:,1,sub) = unifrnd(0,1,[num_of_draws,1]);
        for i=1:num_of_draws
            sv_exp_sixOpt(i,1,sub) = sigma(sub).*((1./cdf_exp_sixOpt(i,1,sub))-1).^(1./beta);
        end
        
        means_exp_sixOpt(sub,1) = mean(sv_exp_sixOpt(:,1,sub));
        stds_exp_sixOpt(sub,1) = std(sv_exp_sixOpt(:,1,sub));

        deviation_mean_s1 = abs((means_exp_sixOpt(sub,1)-mean_sextuplet(1,1,sub))); % compute deviation
        deviation_std_s1 = abs((stds_exp_sixOpt(sub,1)-std_sextuplet(1,1,sub))); % compute deviation
    end
    
    % next draw the small 380 lotteries dist for s2-s6
    deviation_mean = 10000.*ones(6,1); % set arbitrary values
    deviation_std = 10000.*ones(6,1); % set arbitrary values
    
    for j=2:6    
        while deviation_mean(j,1)>threshold_mean_sixOpt(1,j) || deviation_std(j,1)>threshold_std_sixOpt(1,j)
            cdf_exp_sixOpt(:,j,sub) = unifrnd(0,1,[num_of_draws,1]);
            for i=1:length(cdf_exp_sixOpt)
                summationOfSVs=0;
                   for z=1:j-1
                       summationOfSVs = summationOfSVs+(sixOptionsSV(i,z,sub)./sigma(sub)).^beta;
                   end
                sv_exp_sixOpt(i,j,sub) = sigma(sub).*((((1-cdf_exp_sixOpt(i,j,sub)).^(-1./j)-1).^k).*(1+summationOfSVs).^k);
            end
            
            means_exp_sixOpt(sub,j) = mean(sv_exp_sixOpt(:,j,sub));
            stds_exp_sixOpt(sub,j) = std(sv_exp_sixOpt(:,j,sub));

            deviation_mean(j,1) = abs((means_exp_sixOpt(sub,j)-mean_sextuplet(1,j,sub))); % compute deviation
            deviation_std(j,1) = abs((stds_exp_sixOpt(sub,j)-std_sextuplet(1,j,sub))); % compute deviation
        end
    
    end

% export for excel
% SixOpt_exp_table = array2table(sv_exp_sixOpt(:,:,sub),'VariableNames',{'s1','s2','s3','s4','s5','s6'});
% writetable(SixOpt_exp_table,['sv_distributions_sextuplet_sets_sub' num2str(sub) '.csv']);

end

