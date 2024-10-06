%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main code to run through the second stage of the experiment.
% It receives the subjects' choice data from the dirst stage, 
% then estimate the subject-specific utility parameters based on their BDMs.
% Then, it generates the lotteries for the second stage: binary and six-options sets, 
% for both the UNIFORM and PARETO distributions. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% General parameters
num_trials = 320;
max_X_amount = 60; % we don't want to exceed this amount in payments

%% Read first-stage results
% change path according to the saving location of first-stage results file
[subjects_data,text,allraw] = xlsread('sample_first_stage_results.xls');
SID = subjects_data(:,1);
subjects_alpha = subjects_data(:,2);

%% create lotteries for the uniform distribution with binary sets
[binary_choice_DS_uniform, meanSV, max_SV_amount, binary_joint_cool_down_trials] = drawLotteriesSecondStage_binaryUNIFORM(subjects_alpha, num_trials, max_X_amount);

%% create lotteries for the uniform distribution with six-options sets
[six_choice_DS_uniform, six_joint_cool_down_trials] = drawLotteriesSecondStage_sixUNIFORM(subjects_alpha,num_trials,max_X_amount);

%% generate Pareto distributions
[sv_binary_pareto_exp,sv_six_pareto_exp] = pareto_dist_generator_subject_specific_newVersion(subjects_alpha,meanSV(:,1), max_SV_amount);

%% create lotteries for the pareto distribution with binary sets
[binary_choice_DS_pareto] = drawLotteriesSecondStage_binaryPARETO(subjects_alpha,sv_binary_pareto_exp,num_trials,max_X_amount,max_SV_amount, binary_joint_cool_down_trials);

%% create lotteries for the pareto distribution with six-options sets
[six_choice_DS_pareto] = drawLotteriesSecondStage_sixPARETO(subjects_alpha,sv_six_pareto_exp,num_trials,max_X_amount,max_SV_amount,six_joint_cool_down_trials);

%% create subject-specific csv files
for sub=1:length(subjects_alpha)
    % unify cells
    subData = [binary_choice_DS_uniform(:,:,sub) six_choice_DS_uniform(:,:,sub) binary_choice_DS_pareto(:,:,sub) six_choice_DS_pareto(:,:,sub)];  
    fileName = ['sub ' num2str(SID(sub)) '.csv'];
    writecell(subData,fileName);
end



