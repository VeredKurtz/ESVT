function [six_choice_DS] = drawLotteriesSecondStage_sixPARETO(subjects_alpha,sv_six_exp_sixOpt,num_trials,max_X_amount,max_SV_value,joint_cool_down_trials)

%% This function draws the lotteries for the second stage of the ESVT experiment
% We draw 320 pairs of lotteries for binary choice sets.
% For the procedure we use subjects' estimated parameters from the first
% stage of the experiment.
% THIS FUNCTION DRAWS THE SIX-LOTTERIES CHOICE SETS FOR THE *PARETO* DIST CASE 

choise_set_size = 6;
num_joint_trials = 20;
num_of_subjects = length(subjects_alpha);
six_choice_DS = cell(num_trials+1,choise_set_size*4+1,num_of_subjects);


%% SIX OPTIONS CHOICE SET

% figure()
 for s=1:num_of_subjects
    alpha = subjects_alpha(s);
    
    % read PARETO SV files (that were created with the "Pareto generator"
    % function)
%     fileName = ['Output\sv_distributions_sextuplet_sets_sub' num2str(s) '.csv'];
%     distributions=readcell(fileName);
    
    sv(:,1)=sv_six_exp_sixOpt(:,1,s);
    sv(:,2)=sv_six_exp_sixOpt(:,2,s);
    sv(:,3)=sv_six_exp_sixOpt(:,3,s);
    sv(:,4)=sv_six_exp_sixOpt(:,4,s);
    sv(:,5)=sv_six_exp_sixOpt(:,5,s);
    sv(:,6)=sv_six_exp_sixOpt(:,6,s);
    
    % trancuation - find all loterries above max SV value
    for i=1:6
        truncation_locations= find(sv(:,i)>max_SV_value(s));
        sv(truncation_locations,:) = [];
        clear truncation_locations
    end
    
    % we are left with more than 320 SVs, so now pick 320 values at random
    locations = randi(num_trials,num_trials,1);
    for i=1:6
        six_pareto_SV(:,i) = sv(locations,i);
    end
    
    mean_SV(s) = mean(six_pareto_SV(:,1));
%     since we are solving a root when solving for x2, and since x2 cannot
%     be negative, we must define what is the mximal possible x1 we can draw.
%     This is defined by the intersection with the x-axis and the 
%     expression 2SV - X1^alpha.
%     The intersection is at 2SV^(1/alpha).
%     Note that for alpha=1, this converges to 2SV.
%     Plus - we don't want the numbers to be too high, so we need to
%     trancuate our draws up to 60 $.
%     So the max. possible draw of x1 is a minium function: min{2SV^(1/alpha),60}

    max_possible_x1(:,1:choise_set_size) = floor((2.*six_pareto_SV(:,1:choise_set_size)).^(1./alpha));
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
         six_pareto_x2_lottery1(l,1)=100;
         while six_pareto_x2_lottery1(l,1)>max_X_amount
            six_pareto_x1_lottery1(l,1) = round(unifrnd(0,maxX1(l,1)),1); 
            six_pareto_x2_lottery1(l,1) = round(((2.*six_pareto_SV(l,1) - six_pareto_x1_lottery1(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery1 = max(six_pareto_x1_lottery1(l,1),six_pareto_x2_lottery1(l,1));
         min_lottery1 = min(six_pareto_x1_lottery1(l,1),six_pareto_x2_lottery1(l,1));
         six_pareto_x1_lottery1(l,1)= max_lottery1; six_pareto_x2_lottery1(l,1)= min_lottery1;
         
         % LOTTERY II
         six_pareto_x2_lottery2(l,1)=100;
         while six_pareto_x2_lottery2(l,1)>max_X_amount
            six_pareto_x1_lottery2(l,1) = round(unifrnd(0,maxX1(l,2)),1); 
            six_pareto_x2_lottery2(l,1) = round(((2.*six_pareto_SV(l,2) - six_pareto_x1_lottery2(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery2 = max(six_pareto_x1_lottery2(l,1),six_pareto_x2_lottery2(l,1));
         min_lottery2 = min(six_pareto_x1_lottery2(l,1),six_pareto_x2_lottery2(l,1));
         six_pareto_x1_lottery2(l,1)= max_lottery2; six_pareto_x2_lottery2(l,1)= min_lottery2;
         
         % LOTTERY III        
         six_pareto_x2_lottery3(l,1)=100;
         while six_pareto_x2_lottery3(l,1)>max_X_amount
            six_pareto_x1_lottery3(l,1) = round(unifrnd(0,maxX1(l,3)),1); 
            six_pareto_x2_lottery3(l,1) = round(((2.*six_pareto_SV(l,3) - six_pareto_x1_lottery3(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery3 = max(six_pareto_x1_lottery3(l,1),six_pareto_x2_lottery3(l,1));
         min_lottery3 = min(six_pareto_x1_lottery3(l,1),six_pareto_x2_lottery3(l,1));
         six_pareto_x1_lottery3(l,1)= max_lottery3; six_pareto_x2_lottery3(l,1)= min_lottery3;
         
         % LOTTERY IV        
         six_pareto_x2_lottery4(l,1)=100;
         while six_pareto_x2_lottery4(l,1)>max_X_amount
            six_pareto_x1_lottery4(l,1) = round(unifrnd(0,maxX1(l,4)),1); 
            six_pareto_x2_lottery4(l,1) = round(((2.*six_pareto_SV(l,4) - six_pareto_x1_lottery4(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery4 = max(six_pareto_x1_lottery4(l,1),six_pareto_x2_lottery4(l,1));
         min_lottery4 = min(six_pareto_x1_lottery4(l,1),six_pareto_x2_lottery4(l,1));
         six_pareto_x1_lottery4(l,1)= max_lottery4; six_pareto_x2_lottery4(l,1)= min_lottery4;
         
         % LOTTERY V        
         six_pareto_x2_lottery5(l,1)=100;
         while six_pareto_x2_lottery5(l,1)>max_X_amount
            six_pareto_x1_lottery5(l,1) = round(unifrnd(0,maxX1(l,5)),1); 
            six_pareto_x2_lottery5(l,1) = round(((2.*six_pareto_SV(l,5) - six_pareto_x1_lottery5(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery5 = max(six_pareto_x1_lottery5(l,1),six_pareto_x2_lottery5(l,1));
         min_lottery5 = min(six_pareto_x1_lottery5(l,1),six_pareto_x2_lottery5(l,1));
         six_pareto_x1_lottery5(l,1)= max_lottery5; six_pareto_x2_lottery5(l,1)= min_lottery5;
         
         % LOTTERY VI        
         six_pareto_x2_lottery6(l,1)=100;
         while six_pareto_x2_lottery6(l,1)>max_X_amount
            six_pareto_x1_lottery6(l,1) = round(unifrnd(0,maxX1(l,6)),1); 
            six_pareto_x2_lottery6(l,1) = round(((2.*six_pareto_SV(l,6) - six_pareto_x1_lottery6(l,1).^alpha)).^(1/alpha),1);
         end
         % x1 is the high amount
         max_lottery6 = max(six_pareto_x1_lottery6(l,1),six_pareto_x2_lottery6(l,1));
         min_lottery6 = min(six_pareto_x1_lottery6(l,1),six_pareto_x2_lottery6(l,1));
         six_pareto_x1_lottery6(l,1)= max_lottery6; six_pareto_x2_lottery6(l,1)= min_lottery6;
         
             
        % DOMINATION
        % LOTTERY I
        if (six_pareto_x1_lottery1(l,1)<six_pareto_x1_lottery2(l,1)) && (six_pareto_x2_lottery1(l,1)<six_pareto_x2_lottery2(l,1)) 
            fosd_lottery1(l,1) = 1;
        elseif (six_pareto_x1_lottery1(l,1)<six_pareto_x1_lottery3(l,1)) && (six_pareto_x2_lottery1(l,1)<six_pareto_x2_lottery3(l,1)) 
            fosd_lottery1(l,1) = 1;
        elseif (six_pareto_x1_lottery1(l,1)<six_pareto_x1_lottery4(l,1)) && (six_pareto_x2_lottery1(l,1)<six_pareto_x2_lottery4(l,1)) 
            fosd_lottery1(l,1) = 1;
        elseif (six_pareto_x1_lottery1(l,1)<six_pareto_x1_lottery5(l,1)) && (six_pareto_x2_lottery1(l,1)<six_pareto_x2_lottery5(l,1)) 
            fosd_lottery1(l,1) = 1;
        elseif (six_pareto_x1_lottery1(l,1)<six_pareto_x1_lottery6(l,1)) && (six_pareto_x2_lottery1(l,1)<six_pareto_x2_lottery6(l,1)) 
            fosd_lottery1(l,1) = 1;
        else
            fosd_lottery1(l,1) = 0;
        end
        
        % LOTTERY II
        if (six_pareto_x1_lottery2(l,1)<six_pareto_x1_lottery1(l,1)) && (six_pareto_x2_lottery2(l,1)<six_pareto_x2_lottery1(l,1)) 
            fosd_lottery2(l,1) = 1;
        elseif (six_pareto_x1_lottery2(l,1)<six_pareto_x1_lottery3(l,1)) && (six_pareto_x2_lottery2(l,1)<six_pareto_x2_lottery3(l,1)) 
            fosd_lottery2(l,1) = 1;
        elseif (six_pareto_x1_lottery2(l,1)<six_pareto_x1_lottery4(l,1)) && (six_pareto_x2_lottery2(l,1)<six_pareto_x2_lottery4(l,1)) 
            fosd_lottery2(l,1) = 1;
        elseif (six_pareto_x1_lottery2(l,1)<six_pareto_x1_lottery5(l,1)) && (six_pareto_x2_lottery2(l,1)<six_pareto_x2_lottery5(l,1)) 
            fosd_lottery2(l,1) = 1;
        elseif (six_pareto_x1_lottery2(l,1)<six_pareto_x1_lottery6(l,1)) && (six_pareto_x2_lottery2(l,1)<six_pareto_x2_lottery6(l,1)) 
            fosd_lottery2(l,1) = 1;
        else
            fosd_lottery2(l,1) = 0;
        end
        
        % LOTTERY III
        if (six_pareto_x1_lottery3(l,1)<six_pareto_x1_lottery1(l,1)) && (six_pareto_x2_lottery3(l,1)<six_pareto_x2_lottery1(l,1)) 
            fosd_lottery3(l,1) = 1;
        elseif (six_pareto_x1_lottery3(l,1)<six_pareto_x1_lottery2(l,1)) && (six_pareto_x2_lottery3(l,1)<six_pareto_x2_lottery2(l,1)) 
            fosd_lottery3(l,1) = 1;
        elseif (six_pareto_x1_lottery3(l,1)<six_pareto_x1_lottery4(l,1)) && (six_pareto_x2_lottery3(l,1)<six_pareto_x2_lottery4(l,1)) 
            fosd_lottery3(l,1) = 1;
        elseif (six_pareto_x1_lottery3(l,1)<six_pareto_x1_lottery5(l,1)) && (six_pareto_x2_lottery3(l,1)<six_pareto_x2_lottery5(l,1)) 
            fosd_lottery3(l,1) = 1;
        elseif (six_pareto_x1_lottery3(l,1)<six_pareto_x1_lottery6(l,1)) && (six_pareto_x2_lottery3(l,1)<six_pareto_x2_lottery6(l,1)) 
            fosd_lottery3(l,1) = 1;
        else
            fosd_lottery3(l,1) = 0;
        end
        
        % LOTTERY IV
        if (six_pareto_x1_lottery4(l,1)<six_pareto_x1_lottery1(l,1)) && (six_pareto_x2_lottery4(l,1)<six_pareto_x2_lottery1(l,1)) 
            fosd_lottery4(l,1) = 1;
        elseif (six_pareto_x1_lottery4(l,1)<six_pareto_x1_lottery2(l,1)) && (six_pareto_x2_lottery4(l,1)<six_pareto_x2_lottery2(l,1)) 
            fosd_lottery4(l,1) = 1;
        elseif (six_pareto_x1_lottery4(l,1)<six_pareto_x1_lottery3(l,1)) && (six_pareto_x2_lottery4(l,1)<six_pareto_x2_lottery3(l,1)) 
            fosd_lottery4(l,1) = 1;
        elseif (six_pareto_x1_lottery4(l,1)<six_pareto_x1_lottery5(l,1)) && (six_pareto_x2_lottery4(l,1)<six_pareto_x2_lottery5(l,1)) 
            fosd_lottery4(l,1) = 1;
        elseif (six_pareto_x1_lottery4(l,1)<six_pareto_x1_lottery6(l,1)) && (six_pareto_x2_lottery4(l,1)<six_pareto_x2_lottery6(l,1)) 
            fosd_lottery4(l,1) = 1;
        else
            fosd_lottery4(l,1) = 0;
        end
        
        % LOTTERY V
        if (six_pareto_x1_lottery5(l,1)<six_pareto_x1_lottery1(l,1)) && (six_pareto_x2_lottery5(l,1)<six_pareto_x2_lottery1(l,1)) 
            fosd_lottery5(l,1) = 1;
        elseif (six_pareto_x1_lottery5(l,1)<six_pareto_x1_lottery2(l,1)) && (six_pareto_x2_lottery5(l,1)<six_pareto_x2_lottery2(l,1)) 
            fosd_lottery5(l,1) = 1;
        elseif (six_pareto_x1_lottery5(l,1)<six_pareto_x1_lottery3(l,1)) && (six_pareto_x2_lottery5(l,1)<six_pareto_x2_lottery3(l,1)) 
            fosd_lottery5(l,1) = 1;
        elseif (six_pareto_x1_lottery5(l,1)<six_pareto_x1_lottery4(l,1)) && (six_pareto_x2_lottery5(l,1)<six_pareto_x2_lottery4(l,1)) 
            fosd_lottery5(l,1) = 1;
        elseif (six_pareto_x1_lottery5(l,1)<six_pareto_x1_lottery6(l,1)) && (six_pareto_x2_lottery5(l,1)<six_pareto_x2_lottery6(l,1)) 
            fosd_lottery5(l,1) = 1;
        else
            fosd_lottery5(l,1) = 0;
        end
        
        % LOTTERY VI
        if (six_pareto_x1_lottery6(l,1)<six_pareto_x1_lottery1(l,1)) && (six_pareto_x2_lottery6(l,1)<six_pareto_x2_lottery1(l,1)) 
            fosd_lottery6(l,1) = 1;
        elseif (six_pareto_x1_lottery6(l,1)<six_pareto_x1_lottery2(l,1)) && (six_pareto_x2_lottery6(l,1)<six_pareto_x2_lottery2(l,1)) 
            fosd_lottery6(l,1) = 1;
        elseif (six_pareto_x1_lottery6(l,1)<six_pareto_x1_lottery3(l,1)) && (six_pareto_x2_lottery6(l,1)<six_pareto_x2_lottery3(l,1)) 
            fosd_lottery6(l,1) = 1;
        elseif (six_pareto_x1_lottery6(l,1)<six_pareto_x1_lottery4(l,1)) && (six_pareto_x2_lottery6(l,1)<six_pareto_x2_lottery4(l,1)) 
            fosd_lottery6(l,1) = 1;
        elseif (six_pareto_x1_lottery6(l,1)<six_pareto_x1_lottery5(l,1)) && (six_pareto_x2_lottery6(l,1)<six_pareto_x2_lottery5(l,1)) 
            fosd_lottery6(l,1) = 1;
        else
            fosd_lottery6(l,1) = 0;
        end
    end
    
    % sanity check - calculate lotteries' SVs
    six_pareto_sv_lottery1 = ((0.5*six_pareto_x1_lottery1.^alpha + 0.5*six_pareto_x2_lottery1.^alpha));
    six_pareto_sv_lottery2 = ((0.5*six_pareto_x1_lottery2.^alpha + 0.5*six_pareto_x2_lottery2.^alpha));
    six_pareto_sv_lottery3 = ((0.5*six_pareto_x1_lottery3.^alpha + 0.5*six_pareto_x2_lottery3.^alpha));
    six_pareto_sv_lottery4 = ((0.5*six_pareto_x1_lottery4.^alpha + 0.5*six_pareto_x2_lottery4.^alpha));    
    six_pareto_sv_lottery5 = ((0.5*six_pareto_x1_lottery5.^alpha + 0.5*six_pareto_x2_lottery5.^alpha));
    six_pareto_sv_lottery6 = ((0.5*six_pareto_x1_lottery6.^alpha + 0.5*six_pareto_x2_lottery6.^alpha));

    six_pareto_comparison_lotter1 = six_pareto_sv_lottery1-six_pareto_SV(:,1);
    six_pareto_comparison_lotter2 = six_pareto_sv_lottery2-six_pareto_SV(:,2);
    six_pareto_comparison_lotter3 = six_pareto_sv_lottery3-six_pareto_SV(:,3);
    six_pareto_comparison_lotter4 = six_pareto_sv_lottery4-six_pareto_SV(:,4);
    six_pareto_comparison_lotter5 = six_pareto_sv_lottery5-six_pareto_SV(:,5);
    six_pareto_comparison_lotter6 = six_pareto_sv_lottery6-six_pareto_SV(:,6);
    
    six_choice_DS{1,1,s} = 'SV lottery 1 six pareto'; six_choice_DS{1,2,s} = 'SV lottery 2 six pareto'; six_choice_DS{1,3,s} = 'SV lottery 3 six pareto';
    six_choice_DS{1,4,s} = 'SV lottery 4 six pareto'; six_choice_DS{1,5,s} = 'SV lottery 5 six pareto'; six_choice_DS{1,6,s} = 'SV lottery 6 six pareto';
    six_choice_DS{1,7,s} = 'X1 lottery 1 six pareto'; six_choice_DS{1,8,s} = 'X2 lottery 1 six pareto'; six_choice_DS{1,9,s} = 'X1 lottery 2 six pareto';
    six_choice_DS{1,10,s} = 'X2 lottery 2 six pareto'; six_choice_DS{1,11,s} = 'X1 lottery 3 six pareto'; six_choice_DS{1,12,s} = 'X2 lottery 3 six pareto';
    six_choice_DS{1,13,s} = 'X1 lottery 4 six pareto'; six_choice_DS{1,14,s} = 'X2 lottery 4 six pareto'; six_choice_DS{1,15,s} = 'X1 lottery 5 six pareto';
    six_choice_DS{1,16,s} = 'X2 lottery 5 six pareto'; six_choice_DS{1,17,s} = 'X1 lottery 6 six pareto'; six_choice_DS{1,18,s} = 'X2 lottery 6 six pareto';
    six_choice_DS{1,19,s} = 'FOSD lottery 1 six pareto'; six_choice_DS{1,20,s} = 'FOSD lottery 2 six pareto'; six_choice_DS{1,21,s} = 'FOSD lottery 3 six pareto';
    six_choice_DS{1,22,s} = 'FOSD lottery 4 six pareto'; six_choice_DS{1,23,s} = 'FOSD lottery 5 six pareto'; six_choice_DS{1,24,s} = 'FOSD lottery 6 six pareto';
    six_choice_DS{1,25,s} = 'Cool down trial six pareto';
  
    six_choice_DS(2:num_trials+1,1:6,s) = num2cell(six_pareto_SV(:,1:6));
    six_choice_DS(2:num_trials+1,7,s) = num2cell(six_pareto_x1_lottery1(:,1));
    six_choice_DS(2:num_trials+1,8,s) = num2cell(six_pareto_x2_lottery1(:,1));
    six_choice_DS(2:num_trials+1,9,s) = num2cell(six_pareto_x1_lottery2(:,1));
    six_choice_DS(2:num_trials+1,10,s) = num2cell(six_pareto_x2_lottery2(:,1));
    six_choice_DS(2:num_trials+1,11,s) = num2cell(six_pareto_x1_lottery3(:,1));
    six_choice_DS(2:num_trials+1,12,s) = num2cell(six_pareto_x2_lottery3(:,1));
    six_choice_DS(2:num_trials+1,13,s) = num2cell(six_pareto_x1_lottery4(:,1));
    six_choice_DS(2:num_trials+1,14,s) = num2cell(six_pareto_x2_lottery4(:,1));
    six_choice_DS(2:num_trials+1,15,s) = num2cell(six_pareto_x1_lottery5(:,1));
    six_choice_DS(2:num_trials+1,16,s) = num2cell(six_pareto_x2_lottery5(:,1));
    six_choice_DS(2:num_trials+1,17,s) = num2cell(six_pareto_x1_lottery6(:,1));
    six_choice_DS(2:num_trials+1,18,s) = num2cell(six_pareto_x2_lottery6(:,1));
    six_choice_DS(2:num_trials+1,19,s) = num2cell(fosd_lottery1(:,1));
    six_choice_DS(2:num_trials+1,20,s) = num2cell(fosd_lottery2(:,1));
    six_choice_DS(2:num_trials+1,21,s) = num2cell(fosd_lottery3(:,1));
    six_choice_DS(2:num_trials+1,22,s) = num2cell(fosd_lottery4(:,1));
    six_choice_DS(2:num_trials+1,23,s) = num2cell(fosd_lottery5(:,1));
    six_choice_DS(2:num_trials+1,24,s) = num2cell(fosd_lottery6(:,1));


    % randomize cool down trials
    ranomize_cool_down = randperm(num_joint_trials);

    % add joint cool down trials
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,1:6,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,1:6,s));
    
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,7,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,7,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,8,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,8,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,9,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,9,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,10,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,10,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,11,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,11,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,12,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,12,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,13,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,13,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,14,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,14,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,15,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,15,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,16,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,16,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,17,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,17,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,18,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,18,s));
    
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,19,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,19,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,20,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,20,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,21,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,21,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,22,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,22,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,23,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,23,s));
    six_choice_DS(num_trials+2:num_trials+num_joint_trials+1,24,s) = num2cell(joint_cool_down_trials(ranomize_cool_down,24,s));
    
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
    clear truncation_locations_s1 truncation_locations_s2 locations
    clear sv
    clear fosd_lottery1 fosd_lottery2 fosd_lottery3 fosd_lottery4 fosd_lottery5 fosd_lottery6
    
 end

