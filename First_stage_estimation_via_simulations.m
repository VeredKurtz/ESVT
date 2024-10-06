clear
clc

% General parameters
num_trials = 320;
max_X_amount = 60; % we don't want to exceed this amount in payments

%% Read first-stage results
% change path according to the saving location of first-stage results file
dataPath = 'C:\Users\KURTZV02\Dropbox\POST\ESVT\code\session 1\';
[data,text,allraw] = xlsread([dataPath 'session1.xls']);
subjects = unique(data(:,1));
subject_num = length(subjects);
inc = 6/1000000;
simulations_increments = 0:inc:10;

for s=1:subject_num 

    subject = subjects(s);
    subject_rows = find(data(:,1)==subject);
    X(:,1:2) = data(subject_rows,5:6);% lottery data
    ce = data(subject_rows,3); % certainty equivalent
    p=0.5; % equal-prob lotteries;
    
    % start simulations
    for simulation=1:length(simulations_increments)
        alpha = simulations_increments(simulation);
        simulated_CE(:,simulation) = (p.*(X(:,1).^alpha) + p.*(X(:,2).^alpha)).^(1/alpha);
        error_simulations(:,simulation) = (ce - simulated_CE(:,simulation)).^2;
        mean_Euc_Dist(1,simulation) = mean(error_simulations(:,simulation));
        rmse(1,simulation) = sqrt(mean_Euc_Dist(1,simulation));
    end
    
    % find minimal Euclidean distance of all simulations
    [M(s,1),I(s,1)] = min(rmse);
    estimated_alpha(s,1) = simulations_increments(I(s,1));
    min_option = min(X(:,1:2),[],2);
%     if s==4 || s==5
%         figure()
%         scatter(min_option, ce);
%         xlabel('minimal option'); ylabel('CE');
%     end
    
end
    
