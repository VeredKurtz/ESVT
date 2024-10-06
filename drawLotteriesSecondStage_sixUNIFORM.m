function  [six_choice_DS, joint_cool_down_trials] = drawLotteriesSecondStage_sixUNIFORM(subjects_alpha,num_trials,max_X_amount)

%% This function draws the lotteries for the second stage of the ESVT experiment
% We draw 320 sets of six lotteries.
% For the procedure we use subjects' estimated parameters from the first
% stage of the experiment.
% THIS FUNCTION DRAWS THE SIX-LOTTERIES CHOICE SETS FOR THE *UNIFORM* DIST CASE 

choise_set_size = 6;
num_joint_trials = 20;
num_of_subjects = length(subjects_alpha);
six_choice_DS = cell(num_trials+1,choise_set_size*4+1,num_of_subjects);

%% SIX OPTIONS CHOICE SET

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
    rep_sv = repelem(sv_array,8);
    for i=1:choise_set_size
        random_order_sv(:,i) = randperm(length(rep_sv));
        six_unif_SV(:,i) = rep_sv(random_order_sv(:,i));
    end
    
    mean_SV(s) = mean(six_unif_SV(:,1));
%     since we are solving a root when solving for x2, and since x2 cannot
%     be negative, we must define what is the mximal possible x1 we can draw.
%     This is defined by the intersection with the x-axis and the 
%     expression 2SV - X1^alpha.
%     The intersection is at 2SV^(1/alpha).
%     Note that for alpha=1, this converges to 2SV.
%     Plus - we don't want the numbers to be too high, so we need to
%     trancuate our draws up to 60 $.
%     So the max. possible draw of x1 is a minium function: min{2SV^(1/alpha),60}

    max_possible_x1(:,1:choise_set_size) = floor((2.*six_unif_SV(:,1:choise_set_size)).^(1./alpha));
    maxDraw = max_X_amount.*ones(num_trials,1);
    for i=1:choise_set_size
        maxX1(1:num_trials,i) = min([max_possible_x1(:,i) maxDraw],[],2);
    end
    
    % draw the lotteries
    for l=1:num_trials
    %    We draw X1 at random, but its value cannot exceed the max value,
    %    defined above, and then solve for x2.    
    %    We use a power utility function: SV = (p*X1^alpha + p*X2^alpha)
    %    so given p=0.5 --> X2 = (2.*SV - X1^alpha)^(1/alpha)
    %    We make sure to not exceed the maximal possible x2 amount
    
         % LOTTERY I
         six_unif_x2_lottery1(l,1)=100;
         while six_unif_x2_lottery1(l,1)>max_X_amount
            six_unif_x1_lottery1(l,1) = round(unifrnd(0,maxX1(l,1)),1); 
            six_unif_x2_lottery1(l,1) = round(((2.*six_unif_SV(l,1) - six_unif_x1_lottery1(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery1 = max(six_unif_x1_lottery1(l,1),six_unif_x2_lottery1(l,1));
         min_lottery1 = min(six_unif_x1_lottery1(l,1),six_unif_x2_lottery1(l,1));
         six_unif_x1_lottery1(l,1)= max_lottery1; six_unif_x2_lottery1(l,1)= min_lottery1;
         
         % LOTTERY II
         six_unif_x2_lottery2(l,1)=100;
         while six_unif_x2_lottery2(l,1)>max_X_amount
            six_unif_x1_lottery2(l,1) = round(unifrnd(0,maxX1(l,2)),1); 
            six_unif_x2_lottery2(l,1) = round(((2.*six_unif_SV(l,2) - six_unif_x1_lottery2(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery2 = max(six_unif_x1_lottery2(l,1),six_unif_x2_lottery2(l,1));
         min_lottery2 = min(six_unif_x1_lottery2(l,1),six_unif_x2_lottery2(l,1));
         six_unif_x1_lottery2(l,1)= max_lottery2; six_unif_x2_lottery2(l,1)= min_lottery2;
         
         % LOTTERY III        
         six_unif_x2_lottery3(l,1)=100;
         while six_unif_x2_lottery3(l,1)>max_X_amount
            six_unif_x1_lottery3(l,1) = round(unifrnd(0,maxX1(l,3)),1); 
            six_unif_x2_lottery3(l,1) = round(((2.*six_unif_SV(l,3) - six_unif_x1_lottery3(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery3 = max(six_unif_x1_lottery3(l,1),six_unif_x2_lottery3(l,1));
         min_lottery3 = min(six_unif_x1_lottery3(l,1),six_unif_x2_lottery3(l,1));
         six_unif_x1_lottery3(l,1)= max_lottery3; six_unif_x2_lottery3(l,1)= min_lottery3;
         
         % LOTTERY IV        
         six_unif_x2_lottery4(l,1)=100;
         while six_unif_x2_lottery4(l,1)>max_X_amount
            six_unif_x1_lottery4(l,1) = round(unifrnd(0,maxX1(l,4)),1); 
            six_unif_x2_lottery4(l,1) = round(((2.*six_unif_SV(l,4) - six_unif_x1_lottery4(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery4 = max(six_unif_x1_lottery4(l,1),six_unif_x2_lottery4(l,1));
         min_lottery4 = min(six_unif_x1_lottery4(l,1),six_unif_x2_lottery4(l,1));
         six_unif_x1_lottery4(l,1)= max_lottery4; six_unif_x2_lottery4(l,1)= min_lottery4;
         
         % LOTTERY V        
         six_unif_x2_lottery5(l,1)=100;
         while six_unif_x2_lottery5(l,1)>max_X_amount
            six_unif_x1_lottery5(l,1) = round(unifrnd(0,maxX1(l,5)),1); 
            six_unif_x2_lottery5(l,1) = round(((2.*six_unif_SV(l,5) - six_unif_x1_lottery5(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery5 = max(six_unif_x1_lottery5(l,1),six_unif_x2_lottery5(l,1));
         min_lottery5 = min(six_unif_x1_lottery5(l,1),six_unif_x2_lottery5(l,1));
         six_unif_x1_lottery5(l,1)= max_lottery5; six_unif_x2_lottery5(l,1)= min_lottery5;
         
         % LOTTERY VI        
         six_unif_x2_lottery6(l,1)=100;
         while six_unif_x2_lottery6(l,1)>max_X_amount
            six_unif_x1_lottery6(l,1) = round(unifrnd(0,maxX1(l,6)),1); 
            six_unif_x2_lottery6(l,1) = round(((2.*six_unif_SV(l,6) - six_unif_x1_lottery6(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery6 = max(six_unif_x1_lottery6(l,1),six_unif_x2_lottery6(l,1));
         min_lottery6 = min(six_unif_x1_lottery6(l,1),six_unif_x2_lottery6(l,1));
         six_unif_x1_lottery6(l,1)= max_lottery6; six_unif_x2_lottery6(l,1)= min_lottery6;
    
        % DOMINATION
        % LOTTERY I
        if (six_unif_x1_lottery1(l,1)<six_unif_x1_lottery2(l,1)) && (six_unif_x2_lottery1(l,1)<six_unif_x2_lottery2(l,1)) 
            fosd_lottery1(l,1) = 1;
        elseif (six_unif_x1_lottery1(l,1)<six_unif_x1_lottery3(l,1)) && (six_unif_x2_lottery1(l,1)<six_unif_x2_lottery3(l,1)) 
            fosd_lottery1(l,1) = 1;
        elseif (six_unif_x1_lottery1(l,1)<six_unif_x1_lottery4(l,1)) && (six_unif_x2_lottery1(l,1)<six_unif_x2_lottery4(l,1)) 
            fosd_lottery1(l,1) = 1;
        elseif (six_unif_x1_lottery1(l,1)<six_unif_x1_lottery5(l,1)) && (six_unif_x2_lottery1(l,1)<six_unif_x2_lottery5(l,1)) 
            fosd_lottery1(l,1) = 1;
        elseif (six_unif_x1_lottery1(l,1)<six_unif_x1_lottery6(l,1)) && (six_unif_x2_lottery1(l,1)<six_unif_x2_lottery6(l,1)) 
            fosd_lottery1(l,1) = 1;
        else
            fosd_lottery1(l,1) = 0;
        end
        
        % LOTTERY II
        if (six_unif_x1_lottery2(l,1)<six_unif_x1_lottery1(l,1)) && (six_unif_x2_lottery2(l,1)<six_unif_x2_lottery1(l,1)) 
            fosd_lottery2(l,1) = 1;
        elseif (six_unif_x1_lottery2(l,1)<six_unif_x1_lottery3(l,1)) && (six_unif_x2_lottery2(l,1)<six_unif_x2_lottery3(l,1)) 
            fosd_lottery2(l,1) = 1;
        elseif (six_unif_x1_lottery2(l,1)<six_unif_x1_lottery4(l,1)) && (six_unif_x2_lottery2(l,1)<six_unif_x2_lottery4(l,1)) 
            fosd_lottery2(l,1) = 1;
        elseif (six_unif_x1_lottery2(l,1)<six_unif_x1_lottery5(l,1)) && (six_unif_x2_lottery2(l,1)<six_unif_x2_lottery5(l,1)) 
            fosd_lottery2(l,1) = 1;
        elseif (six_unif_x1_lottery2(l,1)<six_unif_x1_lottery6(l,1)) && (six_unif_x2_lottery2(l,1)<six_unif_x2_lottery6(l,1)) 
            fosd_lottery2(l,1) = 1;
        else
            fosd_lottery2(l,1) = 0;
        end
        
        % LOTTERY III
        if (six_unif_x1_lottery3(l,1)<six_unif_x1_lottery1(l,1)) && (six_unif_x2_lottery3(l,1)<six_unif_x2_lottery1(l,1)) 
            fosd_lottery3(l,1) = 1;
        elseif (six_unif_x1_lottery3(l,1)<six_unif_x1_lottery2(l,1)) && (six_unif_x2_lottery3(l,1)<six_unif_x2_lottery2(l,1)) 
            fosd_lottery3(l,1) = 1;
        elseif (six_unif_x1_lottery3(l,1)<six_unif_x1_lottery4(l,1)) && (six_unif_x2_lottery3(l,1)<six_unif_x2_lottery4(l,1)) 
            fosd_lottery3(l,1) = 1;
        elseif (six_unif_x1_lottery3(l,1)<six_unif_x1_lottery5(l,1)) && (six_unif_x2_lottery3(l,1)<six_unif_x2_lottery5(l,1)) 
            fosd_lottery3(l,1) = 1;
        elseif (six_unif_x1_lottery3(l,1)<six_unif_x1_lottery6(l,1)) && (six_unif_x2_lottery3(l,1)<six_unif_x2_lottery6(l,1)) 
            fosd_lottery3(l,1) = 1;
        else
            fosd_lottery3(l,1) = 0;
        end
        
        % LOTTERY IV
        if (six_unif_x1_lottery4(l,1)<six_unif_x1_lottery1(l,1)) && (six_unif_x2_lottery4(l,1)<six_unif_x2_lottery1(l,1)) 
            fosd_lottery4(l,1) = 1;
        elseif (six_unif_x1_lottery4(l,1)<six_unif_x1_lottery2(l,1)) && (six_unif_x2_lottery4(l,1)<six_unif_x2_lottery2(l,1)) 
            fosd_lottery4(l,1) = 1;
        elseif (six_unif_x1_lottery4(l,1)<six_unif_x1_lottery3(l,1)) && (six_unif_x2_lottery4(l,1)<six_unif_x2_lottery3(l,1)) 
            fosd_lottery4(l,1) = 1;
        elseif (six_unif_x1_lottery4(l,1)<six_unif_x1_lottery5(l,1)) && (six_unif_x2_lottery4(l,1)<six_unif_x2_lottery5(l,1)) 
            fosd_lottery4(l,1) = 1;
        elseif (six_unif_x1_lottery4(l,1)<six_unif_x1_lottery6(l,1)) && (six_unif_x2_lottery4(l,1)<six_unif_x2_lottery6(l,1)) 
            fosd_lottery4(l,1) = 1;
        else
            fosd_lottery4(l,1) = 0;
        end
        
        % LOTTERY V
        if (six_unif_x1_lottery5(l,1)<six_unif_x1_lottery1(l,1)) && (six_unif_x2_lottery5(l,1)<six_unif_x2_lottery1(l,1)) 
            fosd_lottery5(l,1) = 1;
        elseif (six_unif_x1_lottery5(l,1)<six_unif_x1_lottery2(l,1)) && (six_unif_x2_lottery5(l,1)<six_unif_x2_lottery2(l,1)) 
            fosd_lottery5(l,1) = 1;
        elseif (six_unif_x1_lottery5(l,1)<six_unif_x1_lottery3(l,1)) && (six_unif_x2_lottery5(l,1)<six_unif_x2_lottery3(l,1)) 
            fosd_lottery5(l,1) = 1;
        elseif (six_unif_x1_lottery5(l,1)<six_unif_x1_lottery4(l,1)) && (six_unif_x2_lottery5(l,1)<six_unif_x2_lottery4(l,1)) 
            fosd_lottery5(l,1) = 1;
        elseif (six_unif_x1_lottery5(l,1)<six_unif_x1_lottery6(l,1)) && (six_unif_x2_lottery5(l,1)<six_unif_x2_lottery6(l,1)) 
            fosd_lottery5(l,1) = 1;
        else
            fosd_lottery5(l,1) = 0;
        end
        
        % LOTTERY VI
        if (six_unif_x1_lottery6(l,1)<six_unif_x1_lottery1(l,1)) && (six_unif_x2_lottery6(l,1)<six_unif_x2_lottery1(l,1)) 
            fosd_lottery6(l,1) = 1;
        elseif (six_unif_x1_lottery6(l,1)<six_unif_x1_lottery2(l,1)) && (six_unif_x2_lottery6(l,1)<six_unif_x2_lottery2(l,1)) 
            fosd_lottery6(l,1) = 1;
        elseif (six_unif_x1_lottery6(l,1)<six_unif_x1_lottery3(l,1)) && (six_unif_x2_lottery6(l,1)<six_unif_x2_lottery3(l,1)) 
            fosd_lottery6(l,1) = 1;
        elseif (six_unif_x1_lottery6(l,1)<six_unif_x1_lottery4(l,1)) && (six_unif_x2_lottery6(l,1)<six_unif_x2_lottery4(l,1)) 
            fosd_lottery6(l,1) = 1;
        elseif (six_unif_x1_lottery6(l,1)<six_unif_x1_lottery5(l,1)) && (six_unif_x2_lottery6(l,1)<six_unif_x2_lottery5(l,1)) 
            fosd_lottery6(l,1) = 1;
        else
            fosd_lottery6(l,1) = 0;
        end
    
    end
    

    % sanity check - calculate lotteries' SVs
    six_unif_sv_lottery1 = ((0.5*six_unif_x1_lottery1.^alpha + 0.5*six_unif_x2_lottery1.^alpha));
    six_unif_sv_lottery2 = ((0.5*six_unif_x1_lottery2.^alpha + 0.5*six_unif_x2_lottery2.^alpha));
    six_unif_sv_lottery3 = ((0.5*six_unif_x1_lottery3.^alpha + 0.5*six_unif_x2_lottery3.^alpha));
    six_unif_sv_lottery4 = ((0.5*six_unif_x1_lottery4.^alpha + 0.5*six_unif_x2_lottery4.^alpha));    
    six_unif_sv_lottery5 = ((0.5*six_unif_x1_lottery5.^alpha + 0.5*six_unif_x2_lottery5.^alpha));
    six_unif_sv_lottery6 = ((0.5*six_unif_x1_lottery6.^alpha + 0.5*six_unif_x2_lottery6.^alpha));

    six_unif_comparison_lotter1 = six_unif_sv_lottery1-six_unif_SV(:,1);
    six_unif_comparison_lotter2 = six_unif_sv_lottery2-six_unif_SV(:,2);
    six_unif_comparison_lotter3 = six_unif_sv_lottery3-six_unif_SV(:,3);
    six_unif_comparison_lotter4 = six_unif_sv_lottery4-six_unif_SV(:,4);
    six_unif_comparison_lotter5 = six_unif_sv_lottery5-six_unif_SV(:,5);
    six_unif_comparison_lotter6 = six_unif_sv_lottery6-six_unif_SV(:,6);
    
    six_choice_DS{1,1,s} = 'SV lottery 1 six uni'; six_choice_DS{1,2,s} = 'SV lottery 2 six uni'; six_choice_DS{1,3,s} = 'SV lottery 3 six uni';
    six_choice_DS{1,4,s} = 'SV lottery 4 six uni'; six_choice_DS{1,5,s} = 'SV lottery 5 six uni'; six_choice_DS{1,6,s} = 'SV lottery 6 six uni';
    six_choice_DS{1,7,s} = 'X1 lottery 1 six uni'; six_choice_DS{1,8,s} = 'X2 lottery 1 six uni'; six_choice_DS{1,9,s} = 'X1 lottery 2 six uni';
    six_choice_DS{1,10,s} = 'X2 lottery 2 six uni'; six_choice_DS{1,11,s} = 'X1 lottery 3 six uni'; six_choice_DS{1,12,s} = 'X2 lottery 3 six uni';
    six_choice_DS{1,13,s} = 'X1 lottery 4 six uni'; six_choice_DS{1,14,s} = 'X2 lottery 4 six uni'; six_choice_DS{1,15,s} = 'X1 lottery 5 six uni';
    six_choice_DS{1,16,s} = 'X2 lottery 5 six uni'; six_choice_DS{1,17,s} = 'X1 lottery 6 six uni'; six_choice_DS{1,18,s} = 'X2 lottery 6 six uni';
    six_choice_DS{1,19,s} = 'FOSD lottery 1 six uni'; six_choice_DS{1,20,s} = 'FOSD lottery 2 six uni'; six_choice_DS{1,21,s} = 'FOSD lottery 3 six uni';
    six_choice_DS{1,22,s} = 'FOSD lottery 4 six uni'; six_choice_DS{1,23,s} = 'FOSD lottery 5 six uni'; six_choice_DS{1,24,s} = 'FOSD lottery 6 six uni';
    six_choice_DS{1,25,s} = 'Cool down trial six uniform';
  
    six_choice_DS(2:num_trials+1,1:6,s) = num2cell(six_unif_SV(:,1:6));
    six_choice_DS(2:num_trials+1,7,s) = num2cell(six_unif_x1_lottery1(:,1));
    six_choice_DS(2:num_trials+1,8,s) = num2cell(six_unif_x2_lottery1(:,1));
    six_choice_DS(2:num_trials+1,9,s) = num2cell(six_unif_x1_lottery2(:,1));
    six_choice_DS(2:num_trials+1,10,s) = num2cell(six_unif_x2_lottery2(:,1));
    six_choice_DS(2:num_trials+1,11,s) = num2cell(six_unif_x1_lottery3(:,1));
    six_choice_DS(2:num_trials+1,12,s) = num2cell(six_unif_x2_lottery3(:,1));
    six_choice_DS(2:num_trials+1,13,s) = num2cell(six_unif_x1_lottery4(:,1));
    six_choice_DS(2:num_trials+1,14,s) = num2cell(six_unif_x2_lottery4(:,1));
    six_choice_DS(2:num_trials+1,15,s) = num2cell(six_unif_x1_lottery5(:,1));
    six_choice_DS(2:num_trials+1,16,s) = num2cell(six_unif_x2_lottery5(:,1));
    six_choice_DS(2:num_trials+1,17,s) = num2cell(six_unif_x1_lottery6(:,1));
    six_choice_DS(2:num_trials+1,18,s) = num2cell(six_unif_x2_lottery6(:,1));
    six_choice_DS(2:num_trials+1,19,s) = num2cell(fosd_lottery1(:,1));
    six_choice_DS(2:num_trials+1,20,s) = num2cell(fosd_lottery2(:,1));
    six_choice_DS(2:num_trials+1,21,s) = num2cell(fosd_lottery3(:,1));
    six_choice_DS(2:num_trials+1,22,s) = num2cell(fosd_lottery4(:,1));
    six_choice_DS(2:num_trials+1,23,s) = num2cell(fosd_lottery5(:,1));
    six_choice_DS(2:num_trials+1,24,s) = num2cell(fosd_lottery6(:,1));


    % Draw joint cool down trials for the end of the session.
    % We cast 20 trials along the diagonal of the SV matrix.
    % One SV is always slightly larger, so these are quite difficult
    % trials.
    % Our SV matrix/array is a 40*40 matrix. We will pick all the even SV
    % bins for this exercise (2nd, 4th, 6th bin, and so forth)
    joint_bins_locations = 2:2:length(sv_array);
    joint_bins_SVs = sv_array(joint_bins_locations);
    
    % Create sextuplets of SVs, such that three are taken from the diagonal and the
    % two SVs are taken from below the diagonal and one SV is taken from
    % above the digonal (to create dominance).
    % 1 - the diagonal; 2 - above; 3 - below
    % Based on our first sample, performence level is at around 75% at a
    % 0.15 distance from the diagonal (where distance is normalized to 0-1
    % scale). In the uniform distribution, this is equivalent to 6th bin
    % away from the diagonal (below or above).
    % Out of the 20 bins we examine - 
    % the first 3 bins - 5 lotteries have to be above from the diagonal and
    % another one from above the diagonal
    % the last 3 bins - 5 lotteries have to be above from the below the diagonal  and
    % the another one from the diagonal.
    below_or_under_diagonal_low    = [ones(3,5) 2.*ones(3,1)];
    below_or_under_diagonal_middle = [ones(length(joint_bins_SVs)-6,3) 2.*ones(length(joint_bins_SVs)-6,1) 3.*ones(length(joint_bins_SVs)-6,2)];
    below_or_under_diagonal_up     = [ones(3,1) 3.*ones(3,5)];
    below_or_under_diagonal        = [below_or_under_diagonal_low; below_or_under_diagonal_middle; below_or_under_diagonal_up];

    for i=1:length(joint_bins_SVs)
        for j = 1:choise_set_size
            if below_or_under_diagonal(i,j)==1
                joint_cool_down_SVs_nonRandom(i,j) = joint_bins_SVs(i);
            elseif below_or_under_diagonal(i,j)==2    
                joint_cool_down_SVs_nonRandom(i,j) = joint_bins_SVs(i)+5.*sv_inc;
            else
                joint_cool_down_SVs_nonRandom(i,j) = joint_bins_SVs(i)-5.*sv_inc;
            end
        end

        % randomize the order of SVs
        rand_order_sv = randperm(choise_set_size);
        for j=1:choise_set_size
            joint_cool_down_SVs(i,j) = joint_cool_down_SVs_nonRandom(i,rand_order_sv(j));
        end
    end

    % find maximal x1
    max_possible_x1_cool_down(:,1:1:choise_set_size) = floor((2.*joint_cool_down_SVs(:,1:1:choise_set_size)).^(1./alpha));
    maxDraw_cool_down = max_X_amount.*ones(num_joint_trials,1);
    for i=1:choise_set_size
        maxX1_cool_down(1:num_joint_trials,i) = min([max_possible_x1_cool_down(:,i) maxDraw_cool_down],[],2);
    end

    % find maximal x1
    max_possible_x1_cool_down(:,1:1:choise_set_size) = floor((2.*joint_cool_down_SVs(:,1:1:choise_set_size)).^(1./alpha));
    maxDraw_cool_down = max_X_amount.*ones(num_joint_trials,1);
    for i=1:choise_set_size
        maxX1_cool_down(1:num_joint_trials,i) = min([max_possible_x1_cool_down(:,i) maxDraw_cool_down],[],2);
    end
    
    % draw joint trials
    for l=1:length(joint_cool_down_SVs)
         % LOTTERY I
         six_unif_x2_lottery1_cool_down(l,1)=100;
         while six_unif_x2_lottery1_cool_down(l,1)>max_X_amount
            six_unif_x1_lottery1_cool_down(l,1) = round(unifrnd(0,maxX1_cool_down(l,1)),1); 
            six_unif_x2_lottery1_cool_down(l,1) = round(((2.*joint_cool_down_SVs(l,1) - six_unif_x1_lottery1_cool_down(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery1_cool_down = max(six_unif_x1_lottery1_cool_down(l,1),six_unif_x2_lottery1_cool_down(l,1));
         min_lottery1_cool_down = min(six_unif_x1_lottery1_cool_down(l,1),six_unif_x2_lottery1_cool_down(l,1));
         six_unif_x1_lottery1_cool_down(l,1)= max_lottery1_cool_down; six_unif_x2_lottery1_cool_down(l,1)= min_lottery1_cool_down;
         
         % LOTTERY II
         six_unif_x2_lottery2_cool_down(l,1)=100;
         while six_unif_x2_lottery2_cool_down(l,1)>max_X_amount
            six_unif_x1_lottery2_cool_down(l,1) = round(unifrnd(0,maxX1_cool_down(l,2)),1); 
            six_unif_x2_lottery2_cool_down(l,1) = round(((2.*joint_cool_down_SVs(l,2) - six_unif_x1_lottery2_cool_down(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery2_cool_down = max(six_unif_x1_lottery2_cool_down(l,1),six_unif_x2_lottery2_cool_down(l,1));
         min_lottery2_cool_down = min(six_unif_x1_lottery2_cool_down(l,1),six_unif_x2_lottery2_cool_down(l,1));
         six_unif_x1_lottery2_cool_down(l,1)= max_lottery2_cool_down; six_unif_x2_lottery2_cool_down(l,1)= min_lottery2_cool_down;
         
         % LOTTERY III        
         six_unif_x2_lottery3_cool_down(l,1)=100;
         while six_unif_x2_lottery3_cool_down(l,1)>max_X_amount
            six_unif_x1_lottery3_cool_down(l,1) = round(unifrnd(0,maxX1_cool_down(l,3)),1); 
            six_unif_x2_lottery3_cool_down(l,1) = round(((2.*joint_cool_down_SVs(l,3) - six_unif_x1_lottery3_cool_down(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery3_cool_down = max(six_unif_x1_lottery3_cool_down(l,1),six_unif_x2_lottery3_cool_down(l,1));
         min_lottery3_cool_down = min(six_unif_x1_lottery3_cool_down(l,1),six_unif_x2_lottery3_cool_down(l,1));
         six_unif_x1_lottery3_cool_down(l,1)= max_lottery3_cool_down; six_unif_x2_lottery3_cool_down(l,1)= min_lottery3_cool_down;
         
         % LOTTERY IV        
         six_unif_x2_lottery4_cool_down(l,1)=100;
         while six_unif_x2_lottery4_cool_down(l,1)>max_X_amount
            six_unif_x1_lottery4_cool_down(l,1) = round(unifrnd(0,maxX1_cool_down(l,4)),1); 
            six_unif_x2_lottery4_cool_down(l,1) = round(((2.*joint_cool_down_SVs(l,4) - six_unif_x1_lottery4_cool_down(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery4_cool_down = max(six_unif_x1_lottery4_cool_down(l,1),six_unif_x2_lottery4_cool_down(l,1));
         min_lottery4_cool_down = min(six_unif_x1_lottery4_cool_down(l,1),six_unif_x2_lottery4_cool_down(l,1));
         six_unif_x1_lottery4_cool_down(l,1)= max_lottery4_cool_down; six_unif_x2_lottery4_cool_down(l,1)= min_lottery4_cool_down;
         
         % LOTTERY V        
         six_unif_x2_lottery5_cool_down(l,1)=100;
         while six_unif_x2_lottery5_cool_down(l,1)>max_X_amount
            six_unif_x1_lottery5_cool_down(l,1) = round(unifrnd(0,maxX1_cool_down(l,5)),1); 
            six_unif_x2_lottery5_cool_down(l,1) = round(((2.*joint_cool_down_SVs(l,5) - six_unif_x1_lottery5_cool_down(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery5_cool_down = max(six_unif_x1_lottery5_cool_down(l,1),six_unif_x2_lottery5_cool_down(l,1));
         min_lottery5_cool_down = min(six_unif_x1_lottery5_cool_down(l,1),six_unif_x2_lottery5_cool_down(l,1));
         six_unif_x1_lottery5_cool_down(l,1)= max_lottery5_cool_down; six_unif_x2_lottery5_cool_down(l,1)= min_lottery5_cool_down;
         
         % LOTTERY VI        
         six_unif_x2_lottery6_cool_down(l,1)=100;
         while six_unif_x2_lottery6_cool_down(l,1)>max_X_amount
            six_unif_x1_lottery6_cool_down(l,1) = round(unifrnd(0,maxX1_cool_down(l,6)),1); 
            six_unif_x2_lottery6_cool_down(l,1) = round(((2.*joint_cool_down_SVs(l,6) - six_unif_x1_lottery6_cool_down(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery6_cool_down = max(six_unif_x1_lottery6_cool_down(l,1),six_unif_x2_lottery6_cool_down(l,1));
         min_lottery6_cool_down = min(six_unif_x1_lottery6_cool_down(l,1),six_unif_x2_lottery6_cool_down(l,1));
         six_unif_x1_lottery6_cool_down(l,1)= max_lottery6_cool_down; six_unif_x2_lottery6_cool_down(l,1)= min_lottery6_cool_down;
    

        % DOMINATION
        % LOTTERY I
        if (six_unif_x1_lottery1_cool_down(l,1)<six_unif_x1_lottery2_cool_down(l,1)) && (six_unif_x2_lottery1_cool_down(l,1)<six_unif_x2_lottery2_cool_down(l,1)) 
            fosd_lottery1_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery1_cool_down(l,1)<six_unif_x1_lottery3_cool_down(l,1)) && (six_unif_x2_lottery1_cool_down(l,1)<six_unif_x2_lottery3_cool_down(l,1)) 
            fosd_lottery1_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery1_cool_down(l,1)<six_unif_x1_lottery4_cool_down(l,1)) && (six_unif_x2_lottery1_cool_down(l,1)<six_unif_x2_lottery4_cool_down(l,1)) 
            fosd_lottery1_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery1_cool_down(l,1)<six_unif_x1_lottery5_cool_down(l,1)) && (six_unif_x2_lottery1_cool_down(l,1)<six_unif_x2_lottery5_cool_down(l,1)) 
            fosd_lottery1_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery1_cool_down(l,1)<six_unif_x1_lottery6_cool_down(l,1)) && (six_unif_x2_lottery1_cool_down(l,1)<six_unif_x2_lottery6_cool_down(l,1)) 
            fosd_lottery1_cool_down(l,1) = 1;
        else
            fosd_lottery1_cool_down(l,1) = 0;
        end
        
        % LOTTERY II
        if (six_unif_x1_lottery2_cool_down(l,1)<six_unif_x1_lottery1_cool_down(l,1)) && (six_unif_x2_lottery2_cool_down(l,1)<six_unif_x2_lottery1_cool_down(l,1)) 
            fosd_lottery2_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery2_cool_down(l,1)<six_unif_x1_lottery3_cool_down(l,1)) && (six_unif_x2_lottery2_cool_down(l,1)<six_unif_x2_lottery3_cool_down(l,1)) 
            fosd_lottery2_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery2_cool_down(l,1)<six_unif_x1_lottery4_cool_down(l,1)) && (six_unif_x2_lottery2_cool_down(l,1)<six_unif_x2_lottery4_cool_down(l,1)) 
            fosd_lottery2_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery2_cool_down(l,1)<six_unif_x1_lottery5_cool_down(l,1)) && (six_unif_x2_lottery2_cool_down(l,1)<six_unif_x2_lottery5_cool_down(l,1)) 
            fosd_lottery2_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery2_cool_down(l,1)<six_unif_x1_lottery6_cool_down(l,1)) && (six_unif_x2_lottery2_cool_down(l,1)<six_unif_x2_lottery6_cool_down(l,1)) 
            fosd_lottery2_cool_down(l,1) = 1;
        else
            fosd_lottery2_cool_down(l,1) = 0;
        end
        
        % LOTTERY III
        if (six_unif_x1_lottery3_cool_down(l,1)<six_unif_x1_lottery1_cool_down(l,1)) && (six_unif_x2_lottery3_cool_down(l,1)<six_unif_x2_lottery1_cool_down(l,1)) 
            fosd_lottery3_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery3_cool_down(l,1)<six_unif_x1_lottery2_cool_down(l,1)) && (six_unif_x2_lottery3_cool_down(l,1)<six_unif_x2_lottery2_cool_down(l,1)) 
            fosd_lottery3_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery3_cool_down(l,1)<six_unif_x1_lottery4_cool_down(l,1)) && (six_unif_x2_lottery3_cool_down(l,1)<six_unif_x2_lottery4_cool_down(l,1)) 
            fosd_lottery3_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery3_cool_down(l,1)<six_unif_x1_lottery5_cool_down(l,1)) && (six_unif_x2_lottery3_cool_down(l,1)<six_unif_x2_lottery5_cool_down(l,1)) 
            fosd_lottery3_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery3_cool_down(l,1)<six_unif_x1_lottery6_cool_down(l,1)) && (six_unif_x2_lottery3_cool_down(l,1)<six_unif_x2_lottery6_cool_down(l,1)) 
            fosd_lottery3_cool_down(l,1) = 1;
        else
            fosd_lottery3_cool_down(l,1) = 0;
        end
        
        % LOTTERY IV
        if (six_unif_x1_lottery4_cool_down(l,1)<six_unif_x1_lottery1_cool_down(l,1)) && (six_unif_x2_lottery4_cool_down(l,1)<six_unif_x2_lottery1_cool_down(l,1)) 
            fosd_lottery4_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery4_cool_down(l,1)<six_unif_x1_lottery2_cool_down(l,1)) && (six_unif_x2_lottery4_cool_down(l,1)<six_unif_x2_lottery2_cool_down(l,1)) 
            fosd_lottery4_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery4_cool_down(l,1)<six_unif_x1_lottery3_cool_down(l,1)) && (six_unif_x2_lottery4_cool_down(l,1)<six_unif_x2_lottery3_cool_down(l,1)) 
            fosd_lottery4_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery4_cool_down(l,1)<six_unif_x1_lottery5_cool_down(l,1)) && (six_unif_x2_lottery4_cool_down(l,1)<six_unif_x2_lottery5_cool_down(l,1)) 
            fosd_lottery4_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery4_cool_down(l,1)<six_unif_x1_lottery6_cool_down(l,1)) && (six_unif_x2_lottery4_cool_down(l,1)<six_unif_x2_lottery6_cool_down(l,1)) 
            fosd_lottery4_cool_down(l,1) = 1;
        else
            fosd_lottery4_cool_down(l,1) = 0;
        end
        
        % LOTTERY V
        if (six_unif_x1_lottery5_cool_down(l,1)<six_unif_x1_lottery1_cool_down(l,1)) && (six_unif_x2_lottery5_cool_down(l,1)<six_unif_x2_lottery1_cool_down(l,1)) 
            fosd_lottery5_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery5_cool_down(l,1)<six_unif_x1_lottery2_cool_down(l,1)) && (six_unif_x2_lottery5_cool_down(l,1)<six_unif_x2_lottery2_cool_down(l,1)) 
            fosd_lottery5_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery5_cool_down(l,1)<six_unif_x1_lottery3_cool_down(l,1)) && (six_unif_x2_lottery5_cool_down(l,1)<six_unif_x2_lottery3_cool_down(l,1)) 
            fosd_lottery5_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery5_cool_down(l,1)<six_unif_x1_lottery4_cool_down(l,1)) && (six_unif_x2_lottery5_cool_down(l,1)<six_unif_x2_lottery4_cool_down(l,1)) 
            fosd_lottery5_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery5_cool_down(l,1)<six_unif_x1_lottery6_cool_down(l,1)) && (six_unif_x2_lottery5_cool_down(l,1)<six_unif_x2_lottery6_cool_down(l,1)) 
            fosd_lottery5_cool_down(l,1) = 1;
        else
            fosd_lottery5_cool_down(l,1) = 0;
        end
        
        % LOTTERY VI
        if (six_unif_x1_lottery6_cool_down(l,1)<six_unif_x1_lottery1_cool_down(l,1)) && (six_unif_x2_lottery6_cool_down(l,1)<six_unif_x2_lottery1_cool_down(l,1)) 
            fosd_lottery6_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery6_cool_down(l,1)<six_unif_x1_lottery2_cool_down(l,1)) && (six_unif_x2_lottery6_cool_down(l,1)<six_unif_x2_lottery2_cool_down(l,1)) 
            fosd_lottery6_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery6_cool_down(l,1)<six_unif_x1_lottery3_cool_down(l,1)) && (six_unif_x2_lottery6_cool_down(l,1)<six_unif_x2_lottery3_cool_down(l,1)) 
            fosd_lottery6_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery6_cool_down(l,1)<six_unif_x1_lottery4_cool_down(l,1)) && (six_unif_x2_lottery6_cool_down(l,1)<six_unif_x2_lottery4_cool_down(l,1)) 
            fosd_lottery6_cool_down(l,1) = 1;
        elseif (six_unif_x1_lottery6_cool_down(l,1)<six_unif_x1_lottery5_cool_down(l,1)) && (six_unif_x2_lottery6_cool_down(l,1)<six_unif_x2_lottery5(l,1)) 
            fosd_lottery6_cool_down(l,1) = 1;
        else
            fosd_lottery6_cool_down(l,1) = 0;
        end
    
    end

    % create output matrix
    joint_cool_down_trials(1:num_joint_trials,1:24,s) = [joint_cool_down_SVs six_unif_x1_lottery1_cool_down six_unif_x2_lottery1_cool_down ...
                                                        six_unif_x1_lottery2_cool_down six_unif_x2_lottery2_cool_down ...
                                                        six_unif_x1_lottery3_cool_down six_unif_x2_lottery3_cool_down ...
                                                        six_unif_x1_lottery4_cool_down six_unif_x2_lottery4_cool_down ...
                                                        six_unif_x1_lottery5_cool_down six_unif_x2_lottery5_cool_down ...
                                                        six_unif_x1_lottery6_cool_down six_unif_x2_lottery6_cool_down ...
                                                        fosd_lottery1_cool_down fosd_lottery2_cool_down fosd_lottery3_cool_down ...
                                                        fosd_lottery4_cool_down fosd_lottery5_cool_down fosd_lottery6_cool_down];

    % randomize cool down trials
    ranomize_cool_down = randperm(num_joint_trials);

    % add cool down trials
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,1:6,s) = num2cell(joint_cool_down_SVs(ranomize_cool_down,1:6));
    
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,7,s) = num2cell(six_unif_x1_lottery1_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,8,s) = num2cell(six_unif_x2_lottery1_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,9,s) = num2cell(six_unif_x1_lottery2_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,10,s) = num2cell(six_unif_x2_lottery2_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,11,s) = num2cell(six_unif_x1_lottery3_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,12,s) = num2cell(six_unif_x2_lottery3_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,13,s) = num2cell(six_unif_x1_lottery4_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,14,s) = num2cell(six_unif_x2_lottery4_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,15,s) = num2cell(six_unif_x1_lottery5_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,16,s) = num2cell(six_unif_x2_lottery5_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,17,s) = num2cell(six_unif_x1_lottery6_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,18,s) = num2cell(six_unif_x2_lottery6_cool_down(ranomize_cool_down,1));
    
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,19,s) = num2cell(fosd_lottery1_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,20,s) = num2cell(fosd_lottery2_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,21,s) = num2cell(fosd_lottery3_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,22,s) = num2cell(fosd_lottery4_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,23,s) = num2cell(fosd_lottery5_cool_down(ranomize_cool_down,1));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,24,s) = num2cell(fosd_lottery6_cool_down(ranomize_cool_down,1));
    
    % add a vector indicating whther this is a cool down trial or not
    six_choice_DS(2:num_trials+1,25,s) = num2cell(0);
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,25,s) = num2cell(1);
    
    clear x1_lottery1 x2_lottery1 sv_lottery1 fileName alpha
    clear max_possible_x1 maxX1
    clear x1_lottery2 x2_lottery2 sv_lottery2  
    clear x1_lottery3 x2_lottery3 sv_lottery3  
    clear x1_lottery4 x2_lottery4 sv_lottery4  
    clear x1_lottery5 x2_lottery5 sv_lottery5  
    clear x1_lottery6 x2_lottery6 sv_lottery6  
    clear fosd_lottery1 fosd_lottery2 fosd_lottery3 fosd_lottery4 fosd_lottery5 fosd_lottery6
    
 end

  