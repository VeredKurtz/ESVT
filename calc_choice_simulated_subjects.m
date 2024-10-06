clear
clc

% read file

data = csvread('simulated six pareto alpha 1.1.csv',1,0);
alpha = 0.3;

%% binary uniform
x1_lottery1_uniform_binary = data(:,3);
x2_lottery1_uniform_binary = data(:,4);
x1_lottery2_uniform_binary = data(:,5);
x2_lottery2_uniform_binary = data(:,6);

expected_utility_lottery1_uniform_binary = 0.5*x1_lottery1_uniform_binary.^alpha + 0.5*x2_lottery1_uniform_binary.^alpha;
expected_utility_lottery2_uniform_binary = 0.5*x1_lottery2_uniform_binary.^alpha + 0.5*x2_lottery2_uniform_binary.^alpha;
EU_binary_uniform = [expected_utility_lottery1_uniform_binary expected_utility_lottery2_uniform_binary];
[~,choices_no_noise_binary_uniform] = max(EU_binary_uniform,[],2);
% add noise in 10% of trials
random_locations_binary_uniform = randi(320,32,1);
[~,min_binary_uniform] = min(EU_binary_uniform(random_locations_binary_uniform,:),[],2);
choices_with_noise_binary_uniform = choices_no_noise_binary_uniform;
choices_with_noise_binary_uniform(random_locations_binary_uniform,1) = min_binary_uniform;

%% binary pareto
x1_lottery1_pareto_binary = data(:,34);
x2_lottery1_pareto_binary = data(:,35);
x1_lottery2_pareto_binary = data(:,36);
x2_lottery2_pareto_binary = data(:,37);

expected_utility_lottery1_pareto_binary = 0.5*x1_lottery1_pareto_binary.^alpha + 0.5*x2_lottery1_pareto_binary.^alpha;
expected_utility_lottery2_pareto_binary = 0.5*x1_lottery2_pareto_binary.^alpha + 0.5*x2_lottery2_pareto_binary.^alpha;
EU_binary_pareto = [expected_utility_lottery1_pareto_binary expected_utility_lottery2_pareto_binary];
[~,choices_no_noise_binary_pareto] = max(EU_binary_pareto,[],2);
% add noise in 10% of trials
random_locations_binary_pareto = randi(320,32,1);
[~,min_binary_pareto] = min(EU_binary_pareto(random_locations_binary_pareto,:),[],2);
choices_with_noise_binary_pareto = choices_no_noise_binary_pareto;
choices_with_noise_binary_pareto(random_locations_binary_pareto,1) = min_binary_pareto;

%% six-options uniform
x1_lottery1_uniform_six = data(:,14); x2_lottery1_uniform_six = data(:,15);
x1_lottery2_uniform_six = data(:,16); x2_lottery2_uniform_six = data(:,17);
x1_lottery3_uniform_six = data(:,18); x2_lottery3_uniform_six = data(:,19);
x1_lottery4_uniform_six = data(:,20); x2_lottery4_uniform_six = data(:,21);
x1_lottery5_uniform_six = data(:,22); x2_lottery5_uniform_six = data(:,23);
x1_lottery6_uniform_six = data(:,24); x2_lottery6_uniform_six = data(:,25);

expected_utility_lottery1_uniform_six = 0.5*x1_lottery1_uniform_six.^alpha + 0.5*x2_lottery1_uniform_six.^alpha;
expected_utility_lottery2_uniform_six = 0.5*x1_lottery2_uniform_six.^alpha + 0.5*x2_lottery2_uniform_six.^alpha;
expected_utility_lottery3_uniform_six = 0.5*x1_lottery3_uniform_six.^alpha + 0.5*x2_lottery3_uniform_six.^alpha;
expected_utility_lottery4_uniform_six = 0.5*x1_lottery4_uniform_six.^alpha + 0.5*x2_lottery4_uniform_six.^alpha;
expected_utility_lottery5_uniform_six = 0.5*x1_lottery5_uniform_six.^alpha + 0.5*x2_lottery5_uniform_six.^alpha;
expected_utility_lottery6_uniform_six = 0.5*x1_lottery6_uniform_six.^alpha + 0.5*x2_lottery6_uniform_six.^alpha;
EU_six_uniform = [expected_utility_lottery1_uniform_six expected_utility_lottery2_uniform_six ...
                  expected_utility_lottery3_uniform_six expected_utility_lottery4_uniform_six ... 
                  expected_utility_lottery5_uniform_six expected_utility_lottery6_uniform_six];
[~,choices_no_noise_six_uniform] = max(EU_six_uniform,[],2);
% add noise in 10% of trials
random_locations_six_uniform = randi(320,32,1);
for i = 1:length(random_locations_binary_pareto)
    index = random_locations_six_uniform(i);
    choice_trial = choices_no_noise_six_uniform(index);
    flag=0;
    while flag==0
        noisy_choice(i,1) = randi(6);
        if noisy_choice(i,1)~=choice_trial
            flag=1;
        else
            flage=0;
        end
    end
end
choices_with_noise_six_uniform = choices_no_noise_six_uniform;
choices_with_noise_six_uniform(random_locations_six_uniform,1) = noisy_choice;

%% six-options pareto
x1_lottery1_pareto_six = data(:,7); x2_lottery1_pareto_six = data(:,8);
x1_lottery2_pareto_six = data(:,9); x2_lottery2_pareto_six = data(:,10);
x1_lottery3_pareto_six = data(:,11); x2_lottery3_pareto_six = data(:,12);
x1_lottery4_pareto_six = data(:,13); x2_lottery4_pareto_six = data(:,14);
x1_lottery5_pareto_six = data(:,15); x2_lottery5_pareto_six = data(:,16);
x1_lottery6_pareto_six = data(:,17); x2_lottery6_pareto_six = data(:,18);

expected_utility_lottery1_pareto_six = 0.5*x1_lottery1_pareto_six.^alpha + 0.5*x2_lottery1_pareto_six.^alpha;
expected_utility_lottery2_pareto_six = 0.5*x1_lottery2_pareto_six.^alpha + 0.5*x2_lottery2_pareto_six.^alpha;
expected_utility_lottery3_pareto_six = 0.5*x1_lottery3_pareto_six.^alpha + 0.5*x2_lottery3_pareto_six.^alpha;
expected_utility_lottery4_pareto_six = 0.5*x1_lottery4_pareto_six.^alpha + 0.5*x2_lottery4_pareto_six.^alpha;
expected_utility_lottery5_pareto_six = 0.5*x1_lottery5_pareto_six.^alpha + 0.5*x2_lottery5_pareto_six.^alpha;
expected_utility_lottery6_pareto_six = 0.5*x1_lottery6_pareto_six.^alpha + 0.5*x2_lottery6_pareto_six.^alpha;
EU_six_pareto = [expected_utility_lottery1_pareto_six expected_utility_lottery2_pareto_six ...
                  expected_utility_lottery3_pareto_six expected_utility_lottery4_pareto_six ... 
                  expected_utility_lottery5_pareto_six expected_utility_lottery6_pareto_six];
[~,choices_no_noise_six_pareto] = max(EU_six_pareto,[],2);
% add noise in 10% of trials
random_locations_six_pareto = randi(320,32,1);
for i = 1:length(random_locations_six_pareto)
    index = random_locations_six_pareto(i);
    choice_trial = choices_no_noise_six_pareto(index);
    flag=0;
    while flag==0
        noisy_choice(i,1) = randi(6);
        if noisy_choice(i,1)~=choice_trial
            flag=1;
        else
            flage=0;
        end
    end
end
choices_with_noise_six_pareto = choices_no_noise_six_pareto;
choices_with_noise_six_pareto(random_locations_six_pareto,1) = noisy_choice;

%% export everything but the six pareto sets
choices = cell(321,6);
choices{1,1} = 'binary uniform no noise';
choices{1,2} = 'six options uniform no noise';
choices{1,3} = 'binary pareto no noise';
choices{1,4} = 'binary uniform with noise';
choices{1,5} = 'six options uniform with noise';
choices{1,6} = 'binary pareto with noise';
choices(2:end,1) = num2cell(choices_no_noise_binary_uniform);
choices(2:end,2) = num2cell(choices_no_noise_six_uniform);
choices(2:end,3) = num2cell(choices_no_noise_binary_pareto);
choices(2:end,4) = num2cell(choices_with_noise_binary_uniform);
choices(2:end,5) = num2cell(choices_with_noise_six_uniform);
choices(2:end,6) = num2cell(choices_with_noise_binary_pareto);
xlswrite('choices_simulated_subject_alpha_1.1.xlsx',choices);

%% export the six pareto sets

choices_six_pareto = cell(321,2);
choices_six_pareto{1,2} = 'six options pareto no noise';
choices_six_pareto{1,5} = 'six options pareto with noise';
choices_six_pareto(2:end,2) = num2cell(choices_no_noise_six_pareto);
choices_six_pareto(2:end,3) = num2cell(choices_with_noise_six_pareto);
xlswrite('choices_six_pareto_simulated_subject_alpha_1.1.xlsx',choices_six_pareto);