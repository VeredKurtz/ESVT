function subjectc_alpha = first_stage_estimation_POWER(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% First stage estimation - ESVT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function estimates the alpha parmeter in a power utility function,
% following the first stage in the ESVT experiment.
% We implement an nls estimation method.
% Last update: 30 April 2021

subjects = unique(data(:,1));
subject_num = length(subjects);

% cast 100 different initial points for the estimation algorithm.
% we will test 100 different initial points to verify the algorithm does
% not reach a local minima
initial_points = unifrnd(0,2,100,1); % test rho values between 0 and 2

for s=1:subject_num 

    subject = subjects(s);
    subject_rows = find(data(:,1)==subject);
    X(:,1:2) = data(subject_rows,5:6);% lottery data
    ce = data(subject_rows,3); % certainty equivalent
    p=0.5; % equal-prob lotteries;
    
    modelfun = @(alpha,x)(p*X(:,1).^(alpha(1))+p*X(:,2).^(alpha(1))).^(1./alpha(1));

    opts = statset('MaxIter',1000);
    coefficients = [];
    for j=1:length(initial_points)
        mdl = fitnlm(X,ce,modelfun,initial_points(j),'CoefficientNames',{['alpha s' num2str(subject)]},...
            'Options',opts); 
        RMSE(j) = mdl.RMSE;
        LogLikelihood(j) = mdl.LogLikelihood;
        AdjRsqr(j) = mdl.Rsquared.Adjusted;
        BIC(j) =mdl.ModelCriterion.BIC;
        coefficients(j,1:4) = table2array(mdl.Coefficients);
    end
    % drop any negative solutions
%     negative_alpha = find(coefficients(:,1)<0);
%     coefficients(negative_alpha,:) = [];
%     RMSE(negative_alpha) = [];
%     BIC(negative_alpha) = [];
%     AdjRsqr(negative_alpha) = [];
%     LogLikelihood(negative_alpha) = [];
     % choose optimal parameters
    [minimal_criterion_POWER, index_POWER] = min(RMSE);
    coefficients_allSubs(s,1:4) = coefficients(index_POWER(1),:);
    RMSE_allSubs(s) = RMSE(index_POWER(1));
    BIC_allSubs(s) = BIC(index_POWER(1));
    AdjRsqr_allSubs(s) = AdjRsqr(index_POWER(1));
    LogLikelihood_allSubs(s) = LogLikelihood(index_POWER(1));
    
    clear x1 x2 ce mdl RMSE LogLikelihood AdjRsqr BIC coeeficients
    
end

%% output
coefficients_table = array2table(coefficients_allSubs,'VariableNames',{'Coeff','SE','tStat','pval'});
RMSE_table = array2table(RMSE_allSubs','VariableNames',{'RMSE'});
LogLikelihood_table = array2table(LogLikelihood_allSubs','VariableNames',{'LogLikelihood'});
AdjRsqr_table = array2table(AdjRsqr_allSubs','VariableNames',{'AdjRsqr'});
BIC_table = array2table(BIC_allSubs','VariableNames',{'BIC'});
subjects_table = array2table(subjects,'VariableNames',{'SubNum'});

output = [subjects_table, coefficients_table, RMSE_table, LogLikelihood_table, AdjRsqr_table, BIC_table];
writetable(output,'firstStageESVT_pilot_powerUtility.csv');

subjectc_alpha = coefficients_allSubs(:,1);
