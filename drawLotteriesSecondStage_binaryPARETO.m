function [binary_choice_DS] = drawLotteriesSecondStage_binaryPARETO(subjects_alpha,sv_binary_pareto_exp,num_trials,max_X_amount,max_SV_value,joint_cool_down_trials)

%% This function draws the lotteries for the second stage of the ESVT experiment
% We draw 320 pairs of lotteries for binary choice sets.
% For the procedure we use subjects' estimated parameters from the first
% stage of the experiment.
% THIS FUNCTION DRAWS THE BINARY-LOTTERIES CHOICE SETS FOR THE *PARETO* DIST CASE 

choise_set_size = 2;
num_joint_trials = 20;
num_of_subjects = length(subjects_alpha);
binary_choice_DS = cell(num_trials+1,(choise_set_size*3)+1,num_of_subjects);
max_FOSDviol = 0.3.*num_trials; % we don't want more than 30% of sets to have FOSD violations
FOSD = zeros(num_trials,1);

%% BINATY CHOICE SET
% figure()
 for sub=1:num_of_subjects
    alpha = subjects_alpha(sub); 
    
    s1=sv_binary_pareto_exp(:,1,sub);
    s2=sv_binary_pareto_exp(:,2,sub);
    
    % trancuation - find all loterries above max SV value
    truncation_locations_s1 = find(s1>max_SV_value(sub));
    s1(truncation_locations_s1) = [];
    s2(truncation_locations_s1) = [];
    truncation_locations_s2 = find(s2>max_SV_value(sub));
    s1(truncation_locations_s2) = [];
    s2(truncation_locations_s2) = [];
      
    % we are left with more than 320 SVs, so now pick 320 values at random
    locations = randi(num_trials,num_trials,1);
    binary_pareto_SV(:,1) = s1(locations(:,1));
    binary_pareto_SV(:,2) = s2(locations(:,1));
  
    % sort by SV difference
    SV_diff = abs(binary_pareto_SV(:,1) - binary_pareto_SV(:,2)); 
    [sv_diff_sort, sv_index] = sort(SV_diff,'descend');
    binary_pareto_SV(:,1:choise_set_size) = binary_pareto_SV(sv_index,1:choise_set_size);
    
    %  Since we are solving a root when solving for x2, and since x2 cannot
    %  be negative, we must define what is the mximal possible x1 we can draw.
    %  This is defined by the intersection with the x-axis and the 
    %  expression 2SV - X1^alpha.
    %  The intersection is at 2SV^(1/alpha).
    %  Note that for alpha=1, this converges to 2SV.
    %  Plus - we don't want the numbers to be too high, so we need to
    %  trancuate our draws up to 60 $.
    %  So the max. possible draw of x1 is a minium function: min{2SV^(1/alpha),60}

    max_possible_x1(:,1:choise_set_size) = floor((2.*binary_pareto_SV(:,1:choise_set_size)).^(1./alpha));
    maxDraw = max_X_amount.*ones(num_trials,1);
    for i=1:choise_set_size
        maxX1(1:num_trials,i) = min([max_possible_x1(:,i) maxDraw],[],2);
    end
    
    % draw the lotteries
    counter_FOSDviol = 0;
    counterFinal = 0;
    for l=1:num_trials
        %  We draw X1 at random, but its value cannot exceed the max value,
        %  defined above, and then solve for x2.    
        %  We use a power utility function: SV = (p*X1^alpha + p*X2^alpha)
        %  so given p=0.5 --> X2 = (2.*SV - X1^alpha)^(1/alpha)
        %  We make sure to not exceed the maximal possible x2 amount
         binary_pareto_x2_lottery1(l,1)=100;
         while binary_pareto_x2_lottery1(l,1)>max_X_amount
            binary_pareto_x1_lottery1(l,1) = round(unifrnd(0,maxX1(l,1)),1); 
            binary_pareto_x2_lottery1(l,1) = round(((2.*binary_pareto_SV(l,1) - binary_pareto_x1_lottery1(l,1).^alpha)).^(1/alpha),1);
         end
         binary_pareto_x2_lottery2(l,1)=100;
         while binary_pareto_x2_lottery2(l,1)>max_X_amount
            binary_pareto_x1_lottery2(l,1) = round(unifrnd(0,maxX1(l,2)),1); 
            binary_pareto_x2_lottery2(l,1) = round(((2.*binary_pareto_SV(l,2) - binary_pareto_x1_lottery2(l,1).^alpha)).^(1/alpha),1);
         end
         
        %  Find if casted lotteries violate FOSD 
        %  First find the higher amount (x1 or x2) in each lottery
        max_lottery1 = max(binary_pareto_x1_lottery1(l,1),binary_pareto_x2_lottery1(l,1));
        min_lottery1 = min(binary_pareto_x1_lottery1(l,1),binary_pareto_x2_lottery1(l,1));
        max_lottery2 = max(binary_pareto_x1_lottery2(l,1),binary_pareto_x2_lottery2(l,1));
        min_lottery2 = min(binary_pareto_x1_lottery2(l,1),binary_pareto_x2_lottery2(l,1));
                 
        if max_lottery1>max_lottery2 && min_lottery1>min_lottery2
            FOSD_viol=2;
        elseif max_lottery2>max_lottery1 && min_lottery2>min_lottery1
            FOSD_viol=1;
        else
            FOSD_viol=0;
        end
        
        % Higher amount is x1
        binary_pareto_x1_lottery1(l,1) = max_lottery1;  binary_pareto_x2_lottery1(l,1) = min_lottery1;
        binary_pareto_x1_lottery2(l,1) = max_lottery2;  binary_pareto_x2_lottery2(l,1) = min_lottery2;
                 
        % Allow for only 45% FOSD violations among the 320 binary lotteries         
        if FOSD_viol>=1
            FOSD(l) = FOSD_viol;
            counter_FOSDviol=counter_FOSDviol+1;
            counterFinal = counterFinal+1;
            if counter_FOSDviol>max_FOSDviol
               % We will re-cast the lottery which has a higher SV
                newCounter = 0;
                while newCounter==0
                    % we have exceeded the amount FOSD lotteries we allow
                    % for in the set.
                    % We therefore cast one of the lotteries again, until
                    % there is no FOSD violation
                    % recast lottery 1
                    binary_pareto_x2_lottery1(l,1)=100;
                     while binary_pareto_x2_lottery1(l,1)>max_X_amount
                        binary_pareto_x1_lottery1(l,1) = round(unifrnd(0,maxX1(l,1)),1); 
                        binary_pareto_x2_lottery1(l,1) = round(((2.*binary_pareto_SV(l,1) - binary_pareto_x1_lottery1(l,1).^alpha)).^(1/alpha),1);
                     end
                     % recast lottery 2
                     binary_pareto_x2_lottery2(l,1)=100;
                     while binary_pareto_x2_lottery2(l,1)>max_X_amount
                        binary_pareto_x1_lottery2(l,1) = round(unifrnd(0,maxX1(l,2)),1); 
                        binary_pareto_x2_lottery2(l,1) = round(((2.*binary_pareto_SV(l,2) - binary_pareto_x1_lottery2(l,1).^alpha)).^(1/alpha),1);
                     end

                    % cheack for FOSD violation of the new binary set
                    max_lottery1 = max(binary_pareto_x1_lottery1(l,1),binary_pareto_x2_lottery1(l,1));
                    min_lottery1 = min(binary_pareto_x1_lottery1(l,1),binary_pareto_x2_lottery1(l,1));
                    max_lottery2 = max(binary_pareto_x1_lottery2(l,1),binary_pareto_x2_lottery2(l,1));
                    min_lottery2 = min(binary_pareto_x1_lottery2(l,1),binary_pareto_x2_lottery2(l,1));
                    FOSD_viol_new = ((max_lottery1>max_lottery2 && min_lottery1>min_lottery2) ...
                         || (max_lottery2>max_lottery1 && min_lottery2>min_lottery1));

                    % if there is no violation, then the counter=1
                    if FOSD_viol_new==0
                        % Higher amount is x1
                        binary_pareto_x1_lottery1(l,1) = max_lottery1;  binary_pareto_x2_lottery1(l,1) = min_lottery1;
                        binary_pareto_x1_lottery2(l,1) = max_lottery2;  binary_pareto_x2_lottery2(l,1) = min_lottery2;
                        % keep track we do not exceed max violations
                        % allowed
                        counterFinal = counterFinal-1;
                        FOSD(l) = 0;
                        newCounter=1;
                    end
                end    
            end
            clear high_SV
        end
        clear FOSD_viol FOSD_viol_new
    end
           
    %   sanity check - calculate lotteries' SVs
    binary_pareto_sv_lottery1 = ((0.5*binary_pareto_x1_lottery1.^alpha + 0.5*binary_pareto_x2_lottery1.^alpha));
    binary_pareto_sv_lottery2 = ((0.5*binary_pareto_x1_lottery2.^alpha + 0.5*binary_pareto_x2_lottery2.^alpha));

    binary_pareto_comparison_lotter1 = binary_pareto_sv_lottery1-binary_pareto_SV(:,1);
    binary_pareto_comparison_lotter2 = binary_pareto_sv_lottery2-binary_pareto_SV(:,2);

    binary_choice_DS{1,1,sub} = 'SV lottery 1 binary pareto'; 
    binary_choice_DS{1,2,sub} = 'SV lottery 2 binary pareto';
    binary_choice_DS{1,3,sub} = 'X1 lottery 1 binary pareto';
    binary_choice_DS{1,4,sub} = 'X2 lottery 1 binary pareto';
    binary_choice_DS{1,5,sub} = 'X1 lottery 2 binary pareto';
    binary_choice_DS{1,6,sub} = 'X2 lottery 2 binary pareto';
    binary_choice_DS{1,7,sub} = 'FOSD binary pareto';
    binary_choice_DS{1,8,sub} = 'Cool down trial binary pareto';
    
    % randomize trials
    randomizationOfTrials = randperm(num_trials);
    % save in one matrix
    binary_choice_DS(2:num_trials+1,1:2,sub) = num2cell(binary_pareto_SV(randomizationOfTrials,1:2));
    binary_choice_DS(2:num_trials+1,3,sub) = num2cell(binary_pareto_x1_lottery1(randomizationOfTrials,1));
    binary_choice_DS(2:num_trials+1,4,sub) = num2cell(binary_pareto_x2_lottery1(randomizationOfTrials,1));
    binary_choice_DS(2:num_trials+1,5,sub) = num2cell(binary_pareto_x1_lottery2(randomizationOfTrials,1));
    binary_choice_DS(2:num_trials+1,6,sub) = num2cell(binary_pareto_x2_lottery2(randomizationOfTrials,1));
    binary_choice_DS(2:num_trials+1,7,sub) = num2cell(FOSD(randomizationOfTrials));

    
    % total number of FOSD violations
    sumFOSD(sub,1) = counterFinal;
    % Find means of the distributions for the PARETO case
    meanSV(sub,1:2) = mean(binary_pareto_SV(:,1:2));
       
    % add joint cool down trials
    % randomize cool down trials
    ranomize_cool_down = randperm(length(joint_cool_down_trials));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,1:2,sub) = num2cell(joint_cool_down_trials(ranomize_cool_down,1:2,sub));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,3,sub) = num2cell(joint_cool_down_trials(ranomize_cool_down,3,sub));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,4,sub) = num2cell(joint_cool_down_trials(ranomize_cool_down,4,sub));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,5,sub) = num2cell(joint_cool_down_trials(ranomize_cool_down,5,sub));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,6,sub) = num2cell(joint_cool_down_trials(ranomize_cool_down,6,sub));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,7,sub) = num2cell(joint_cool_down_trials(ranomize_cool_down,7,sub));
    
    % add a vector indicating whther this is a cool down trial or not
    binary_choice_DS(2:num_trials+1,8,sub) = num2cell(0);
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,8,sub) = num2cell(1);

    clear alpha x1_lottery1 x2_lottery1 sv_lottery1 sv_lottery2 fileName
    clear FOSD_viol x1_lottery2 x2_lottery2 max_possible_x1 maxX1 counterFinal counter_FOSDviol  
    clear truncation_locations_s1 truncation_locations_s2 locations
    
 end



