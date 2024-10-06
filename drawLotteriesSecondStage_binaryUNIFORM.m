function [binary_choice_DS, meanSV, max_SV_amount, joint_cool_down_trials]  = drawLotteriesSecondStage_binaryUNIFORM(subjects_alpha,num_trials,max_X_amount)

%% This function draws the lotteries for the second stage of the ESVT experiment
% We draw 320 pairs of lotteries for binary choice sets.
% For the procedure we use subjects' estimated parameters from the first
% stage of the experiment.
% THIS FUNCTION DRAWS THE BINARY-LOTTERIES CHOICE SETS FOR THE *UNIFORM* DIST CASE 

%%
choise_set_size = 2;
num_of_subjects = length(subjects_alpha);
num_joint_trials = 20;
binary_choice_DS = cell(num_trials+num_joint_trials+1,(choise_set_size*3)+2,num_of_subjects);
max_FOSDviol = 0.45.*num_trials; % we don't want more than 45% of sets to have FOSD violations
FOSD = zeros(num_trials,1);

%% BINATY CHOICE SET

% figure()
 for s=1:num_of_subjects
    alpha = subjects_alpha(s); 
    % the maximal possible SV would be a lottery in which both Xs have the
    % maximal possible value (60$)
    max_SV_amount(s) = 0.5.*max_X_amount.^alpha + 0.5.*max_X_amount.^alpha;
    % Create a uniform distribution of SV with 40 different values.
    % Each SV value will be repeated 8 times.
    sv_inc = (max_SV_amount(s)-1)./39;
    sv_array = 1:sv_inc:max_SV_amount(s);
    sv_array = sv_array';
    repetitions = num_trials./length(sv_array);
    rep_sv = repelem(sv_array,repetitions);
    for i=1:choise_set_size
        random_order_sv(:,i) = randperm(length(rep_sv));
        binary_unif_SV(:,i) = rep_sv(random_order_sv(:,i));
    end
    
    % sort by SV difference
    SV_diff = abs(binary_unif_SV(:,1) - binary_unif_SV(:,2)); 
    [sv_diff_sort, sv_index] = sort(SV_diff,'descend');
    binary_unif_SV(:,1:choise_set_size) = binary_unif_SV(sv_index,1:choise_set_size);
    
%     since we are solving a root when solving for x2, and since x2 cannot
%     be negative, we must define what is the mximal possible x1 we can draw.
%     This is defined by the intersection with the x-axis and the 
%     expression 2SV - X1^alpha.
%     The intersection is at 2SV^(1/alpha).
%     Note that for alpha=1, this converges to 2SV.
%     Plus - we don't want the numbers to be too high, so we need to
%     trancuate our draws up to 60 $.
%     So the max. possible draw of x1 is a minium function: min{2SV^(1/alpha),60}

    max_possible_x1(:,1:choise_set_size) = floor((2.*binary_unif_SV(:,1:choise_set_size)).^(1./alpha));
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
         binary_unif_x2_lottery1(l,1)=100;
         while binary_unif_x2_lottery1(l,1)>max_X_amount
            binary_unif_x1_lottery1(l,1) = round(unifrnd(0,maxX1(l,1)),1); 
            binary_unif_x2_lottery1(l,1) = round(((2.*binary_unif_SV(l,1) - binary_unif_x1_lottery1(l,1).^alpha)).^(1/alpha),1);
         end
         binary_unif_x2_lottery2(l,1)=100;
         while binary_unif_x2_lottery2(l,1)>max_X_amount
            binary_unif_x1_lottery2(l,1) = round(unifrnd(0,maxX1(l,2)),1); 
            binary_unif_x2_lottery2(l,1) = round(((2.*binary_unif_SV(l,2) - binary_unif_x1_lottery2(l,1).^alpha)).^(1/alpha),1);
         end
         
        %  Find if casted lotteries violate FOSD 
        %  First find the higher amount (x1 or x2) in each lottery
        max_lottery1 = max(binary_unif_x1_lottery1(l,1),binary_unif_x2_lottery1(l,1));
        min_lottery1 = min(binary_unif_x1_lottery1(l,1),binary_unif_x2_lottery1(l,1));
        max_lottery2 = max(binary_unif_x1_lottery2(l,1),binary_unif_x2_lottery2(l,1));
        min_lottery2 = min(binary_unif_x1_lottery2(l,1),binary_unif_x2_lottery2(l,1));
        if max_lottery1>max_lottery2 && min_lottery1>min_lottery2
            FOSD_viol=2;
        elseif max_lottery2>max_lottery1 && min_lottery2>min_lottery1
            FOSD_viol=1;
        else
            FOSD_viol=0;
        end
        
        % Higher amount is x1
        binary_unif_x1_lottery1(l,1) = max_lottery1;  binary_unif_x2_lottery1(l,1) = min_lottery1;
        binary_unif_x1_lottery2(l,1) = max_lottery2;  binary_unif_x2_lottery2(l,1) = min_lottery2;
                 
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
                    binary_unif_x2_lottery1(l,1)=100;
                     while binary_unif_x2_lottery1(l,1)>max_X_amount
                        binary_unif_x1_lottery1(l,1) = round(unifrnd(0,maxX1(l,1)),1); 
                        binary_unif_x2_lottery1(l,1) = round(((2.*binary_unif_SV(l,1) - binary_unif_x1_lottery1(l,1).^alpha)).^(1/alpha),1);
                     end
                     % recast lottery 2
                     binary_unif_x2_lottery2(l,1)=100;
                     while binary_unif_x2_lottery2(l,1)>max_X_amount
                        binary_unif_x1_lottery2(l,1) = round(unifrnd(0,maxX1(l,2)),1); 
                        binary_unif_x2_lottery2(l,1) = round(((2.*binary_unif_SV(l,2) - binary_unif_x1_lottery2(l,1).^alpha)).^(1/alpha),1);
                     end

                    % check for FOSD violation of the new binary set
                    max_lottery1 = max(binary_unif_x1_lottery1(l,1),binary_unif_x2_lottery1(l,1));
                    min_lottery1 = min(binary_unif_x1_lottery1(l,1),binary_unif_x2_lottery1(l,1));
                    max_lottery2 = max(binary_unif_x1_lottery2(l,1),binary_unif_x2_lottery2(l,1));
                    min_lottery2 = min(binary_unif_x1_lottery2(l,1),binary_unif_x2_lottery2(l,1));
                    FOSD_viol_new = ((max_lottery1>max_lottery2 && min_lottery1>min_lottery2) ...
                         || (max_lottery2>max_lottery1 && min_lottery2>min_lottery1));
                     
                    % if there is no violation, then the counter=1
                    if FOSD_viol_new==0
                        % Higher amount is x1
                        binary_unif_x1_lottery1(l,1) = max_lottery1;  binary_unif_x2_lottery1(l,1) = min_lottery1;
                        binary_unif_x1_lottery2(l,1) = max_lottery2;  binary_unif_x2_lottery2(l,1) = min_lottery2;
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
    binary_unif_sv_lottery1 = ((0.5*binary_unif_x1_lottery1.^alpha + 0.5*binary_unif_x2_lottery1.^alpha));
    binary_unif_sv_lottery2 = ((0.5*binary_unif_x1_lottery2.^alpha + 0.5*binary_unif_x2_lottery2.^alpha));

    binary_unif_comparison_lotter1 = binary_unif_sv_lottery1-binary_unif_SV(:,1);
    binary_unif_comparison_lotter2 = binary_unif_sv_lottery2-binary_unif_SV(:,2);

    binary_choice_DS{1,1,s} = 'SV lottery 1 binary uni'; 
    binary_choice_DS{1,2,s} = 'SV lottery 2 binary uni';
    binary_choice_DS{1,3,s} = 'X1 lottery 1 binary uni';
    binary_choice_DS{1,4,s} = 'X2 lottery 1 binary uni';
    binary_choice_DS{1,5,s} = 'X1 lottery 2 binary uni';
    binary_choice_DS{1,6,s} = 'X2 lottery 2 binary uni';
    binary_choice_DS{1,7,s} = 'FOSD binary uniform';
    binary_choice_DS{1,8,s} = 'Cool down trial binary uniform';
    
    % randomize trials
    randomizationOfTrials = randperm(num_trials);
    % save in one matrix
    binary_choice_DS(2:num_trials+1,1:2,s) = num2cell(binary_unif_SV(randomizationOfTrials,1:2));
    binary_choice_DS(2:num_trials+1,3,s) = num2cell(binary_unif_x1_lottery1(randomizationOfTrials,1));
    binary_choice_DS(2:num_trials+1,4,s) = num2cell(binary_unif_x2_lottery1(randomizationOfTrials,1));
    binary_choice_DS(2:num_trials+1,5,s) = num2cell(binary_unif_x1_lottery2(randomizationOfTrials,1));
    binary_choice_DS(2:num_trials+1,6,s) = num2cell(binary_unif_x2_lottery2(randomizationOfTrials,1));
    binary_choice_DS(2:num_trials+1,7,s) = num2cell(FOSD(randomizationOfTrials));
    
    % total number of FOSD violations
    sumFOSD(s,1) = counterFinal;
    % Find means of the distributions for the PARETO case
    meanSV(s,1:2) = mean(binary_unif_SV(:,1:2));
       
    
    % Draw joint cool down trials for the end of the session.
    % We cast 20 trials along the diagonal of the SV matrix.
    % One SV is always slightly larger, so these are quite difficult
    % trials.
    % Our SV matrix/array is a 40*40 matrix. We will pick all the even SV
    % bins for this exercise (2nd, 4th, 6th bin, and so forth)
    joint_bins_locations = 2:2:length(sv_array);
    joint_bins_SVs = sv_array(joint_bins_locations);

    % Create pairs of SVs, such that one is taken from the diagonal and the
    % second one is either alternating from above the diagonal or below it.
    % Based on our first sample, performence level is at around 75% at a
    % 0.15 distance from the diagonal (where distance is normalized to 0-1
    % scale). In the uniform distribution, this is equivalent to 6th bin
    % away from the diagonal (below or above).
    % Out of the 20 bins we examine - 
    % the first 3 bins have to be above the diagonal.
    % the last 3 bins have to be belpw the diagonal
    % the rest we cast at random.
    below_or_under_diagonal = [1; 1; 1; randi(2,length(joint_bins_SVs)-6,1); 2; 2; 2];
    for i=1:length(joint_bins_SVs)
        joint_cool_down_SVs_nonRandom(i,1) = joint_bins_SVs(i);
        if below_or_under_diagonal(i,1)==1
            joint_cool_down_SVs_nonRandom(i,2) = joint_bins_SVs(i)+5.*sv_inc;
        else
            joint_cool_down_SVs_nonRandom(i,2) = joint_bins_SVs(i)-5.*sv_inc;
        end
        % randomize the order of the two SVs
        rand_order_sv = randperm(2);
        for j=1:2
            joint_cool_down_SVs(i,j) = joint_cool_down_SVs_nonRandom(i,rand_order_sv(j));
        end
    end

    % find maximal x1
    max_possible_x1_cool_down(:,1:1:choise_set_size) = floor((2.*joint_cool_down_SVs(:,1:1:choise_set_size)).^(1./alpha));
    maxDraw_cool_down = max_X_amount.*ones(num_joint_trials,1);
    for i=1:choise_set_size
        maxX1_cool_down(1:num_joint_trials,i) = min([max_possible_x1_cool_down(:,i) maxDraw_cool_down],[],2);
    end
    
    % draw trials
    for j=1:length(joint_cool_down_SVs)
        binary_unif_x2_lottery1_cool_down(j,1)=100;
        while binary_unif_x2_lottery1_cool_down(j,1)>max_X_amount
            binary_unif_x1_lottery1_cool_down(j,1) = round(unifrnd(0,maxX1_cool_down(j,1)),1); 
            binary_unif_x2_lottery1_cool_down(j,1) = round(((2.*joint_cool_down_SVs(j,1) - binary_unif_x1_lottery1_cool_down(j,1).^alpha)).^(1/alpha),1);
        end
        binary_unif_x2_lottery2_cool_down(j,1)=100;
        while binary_unif_x2_lottery2_cool_down(j,1)>max_X_amount
            binary_unif_x1_lottery2_cool_down(j,1) = round(unifrnd(0,maxX1_cool_down(j,2)),1); 
            binary_unif_x2_lottery2_cool_down(j,1) = round(((2.*joint_cool_down_SVs(j,2) - binary_unif_x1_lottery2_cool_down(j,1).^alpha)).^(1/alpha),1);
        end

        % check for FOSD violation of the binary set
        max_lottery1 = max(binary_unif_x1_lottery1_cool_down(j,1),binary_unif_x2_lottery1_cool_down(j,1));
        min_lottery1 = min(binary_unif_x1_lottery1_cool_down(j,1),binary_unif_x2_lottery1_cool_down(j,1));
        max_lottery2 = max(binary_unif_x1_lottery2_cool_down(j,1),binary_unif_x2_lottery2_cool_down(j,1));
        min_lottery2 = min(binary_unif_x1_lottery2_cool_down(j,1),binary_unif_x2_lottery2_cool_down(j,1));
        FOSD_viol = ((max_lottery1>max_lottery2 && min_lottery1>min_lottery2) ...
             || (max_lottery2>max_lottery1 && min_lottery2>min_lottery1));
         
        if FOSD_viol==0
            FOSD_cool_down(j,1) = 0; 
        else
            FOSD_cool_down(j,1) = 1; 
        end

        % Higher amount is x1
        binary_unif_x1_lottery1_cool_down(j,1) = max_lottery1;  binary_unif_x2_lottery1_cool_down(j,1) = min_lottery1;
        binary_unif_x1_lottery2_cool_down(j,1) = max_lottery2;  binary_unif_x2_lottery2_cool_down(j,1) = min_lottery2;    

        % recast trials with FOSD 
        if FOSD_viol>=1 
           % We will re-cast the lottery which has a higher SV
            newCounter = 0;
            while newCounter==0.
                % We cast one of the lotteries again, until
                % there is no FOSD violation
                % recast lottery 1
                binary_unif_x2_lottery1_cool_down(j,1)=100;
                 while binary_unif_x2_lottery1_cool_down(j,1)>max_X_amount
                    binary_unif_x1_lottery1_cool_down(j,1) = round(unifrnd(0,maxX1_cool_down(j,1)),1); 
                    binary_unif_x2_lottery1_cool_down(j,1) = round(((2.*joint_cool_down_SVs(j,1) - binary_unif_x1_lottery1_cool_down(j,1).^alpha)).^(1/alpha),1);
                 end
                 % recast lottery 2
                 binary_unif_x2_lottery2_cool_down(j,1)=100;
                 while binary_unif_x2_lottery2_cool_down(j,1)>max_X_amount
                    binary_unif_x1_lottery2_cool_down(j,1) = round(unifrnd(0,maxX1_cool_down(j,2)),1); 
                    binary_unif_x2_lottery2_cool_down(j,1) = round(((2.*joint_cool_down_SVs(j,2) - binary_unif_x1_lottery2_cool_down(j,1).^alpha)).^(1/alpha),1);
                 end

                % check for FOSD violation of the new binary set
                max_lottery1 = max(binary_unif_x1_lottery1_cool_down(j,1),binary_unif_x2_lottery1_cool_down(j,1));
                min_lottery1 = min(binary_unif_x1_lottery1_cool_down(j,1),binary_unif_x2_lottery1_cool_down(j,1));
                max_lottery2 = max(binary_unif_x1_lottery2_cool_down(j,1),binary_unif_x2_lottery2_cool_down(j,1));
                min_lottery2 = min(binary_unif_x1_lottery2_cool_down(j,1),binary_unif_x2_lottery2_cool_down(j,1));
                FOSD_viol_new = ((max_lottery1>max_lottery2 && min_lottery1>min_lottery2) ...
                     || (max_lottery2>max_lottery1 && min_lottery2>min_lottery1));
                 
                % if there is no violation, then the counter=1
                if FOSD_viol_new==0
                    % Higher amount is x1
                    binary_unif_x1_lottery1_cool_down(j,1) = max_lottery1;  binary_unif_x2_lottery1_cool_down(j,1) = min_lottery1;
                    binary_unif_x1_lottery2_cool_down(j,1) = max_lottery2;  binary_unif_x2_lottery2_cool_down(j,1) = min_lottery2;    
                    FOSD_cool_down(j,1) = 0; 
                    newCounter=1;
                end
            end
            
            clear FOSD_viol FOSD_viol_new
        end
    
    end
    
    % create output matrix
    joint_cool_down_trials(1:num_joint_trials,1:7,s) = [joint_cool_down_SVs binary_unif_x1_lottery1_cool_down binary_unif_x2_lottery1_cool_down binary_unif_x1_lottery2_cool_down binary_unif_x2_lottery2_cool_down FOSD_cool_down];

    % randomize cool down trials
    ranomize_cool_down = randperm(length(joint_cool_down_trials));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,1:2,s) = num2cell(joint_cool_down_SVs(ranomize_cool_down,1:2));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,3,s) = num2cell(binary_unif_x1_lottery1_cool_down(ranomize_cool_down,1));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,4,s) = num2cell(binary_unif_x2_lottery1_cool_down(ranomize_cool_down,1));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,5,s) = num2cell(binary_unif_x1_lottery2_cool_down(ranomize_cool_down,1));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,6,s) = num2cell(binary_unif_x2_lottery2_cool_down(ranomize_cool_down,1));
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,7,s) = num2cell(FOSD_cool_down(ranomize_cool_down,1));
    
    % add a vector indicating whther this is a cool down trial or not
    binary_choice_DS(2:num_trials+1,8,s) = num2cell(0);
    binary_choice_DS(num_trials+2:num_trials+length(joint_cool_down_trials)+1,8,s) = num2cell(1);

    clear alpha x1_lottery1 x2_lottery1 sv_lottery1 sv_lottery2 fileName
    clear FOSD_viol x1_lottery2 x2_lottery2 max_possible_x1 maxX1 counterFinal counter_FOSDviol   
    
 end

