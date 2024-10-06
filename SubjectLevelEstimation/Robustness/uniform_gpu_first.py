import random
import pandas as pd
import torch
import torch.nn.functional as F
import numpy as np
from scipy.optimize import minimize
import pandas as pd 
import numpy as np
import logging

logging.basicConfig(format='%(asctime)s,%(msecs)03d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
    datefmt='%Y-%m-%d:%H:%M:%S',
    level=logging.DEBUG)


def obj_function(parameters, option1_high, option1_low, option2_high, option2_low, actual_chosen_option):


    num_simulations_per_trial = 10000
    num_trials = len(option1_high)
    alpha, M, delta = parameters
    # delta = 0.03

    actual_chosen_option = torch.from_numpy(actual_chosen_option).to('cuda').float()
    option1_high = torch.from_numpy((np.repeat(option1_high, num_simulations_per_trial)).reshape(num_trials, num_simulations_per_trial)).to('cuda').float()
    option1_low = torch.from_numpy((np.repeat(option1_low, num_simulations_per_trial)).reshape(num_trials, num_simulations_per_trial)).to('cuda').float()
    option2_high = torch.from_numpy((np.repeat(option2_high, num_simulations_per_trial)).reshape(num_trials, num_simulations_per_trial)).to('cuda').float()
    option2_low = torch.from_numpy((np.repeat(option2_low, num_simulations_per_trial)).reshape(num_trials, num_simulations_per_trial)).to('cuda').float()
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    N_option1_high = ((option1_high**alpha) / ((option1_high**alpha) + (M**alpha))) 
    N_option1_low = ((option1_low**alpha) / ((option1_low**alpha) + (M**alpha))) 
    N_option2_high = ((option2_high**alpha) / ((option2_high**alpha) + (M**alpha))) 
    N_option2_low = ((option2_low**alpha) / ((option2_low**alpha) + (M**alpha))) 

    N_option1 = (0.5 * N_option1_high) + (0.5 * N_option1_low)
    N_option2 = (0.5 * N_option2_high) + (0.5 * N_option2_low)

    mean_tensor = torch.full((num_trials, num_simulations_per_trial), 0.0, device=device)
    std_dev_tensor = torch.full((num_trials, num_simulations_per_trial), delta, device=device)

    N_option1 = N_option1 + torch.normal(mean=mean_tensor, std=std_dev_tensor)
    N_option2 = N_option2 + torch.normal(mean=mean_tensor, std=std_dev_tensor)

    stacked_options = torch.stack((N_option1, N_option2), dim=2)
    chosen_option2 = torch.argmax(stacked_options, dim=2).float()
    chosen_option1 = 1 - chosen_option2

    choice_prob_option1 = torch.mean(chosen_option1, dim=1)
    choice_prob_option2 = torch.mean(chosen_option2, dim=1)

    small_value = 1e-5
    likelihood = (actual_chosen_option * choice_prob_option2) + ((1 - actual_chosen_option) * choice_prob_option1)
    likelihood[likelihood == 0] = small_value

    log_likelihood = torch.sum(torch.log(likelihood))

    return -log_likelihood.item()  

all_data = pd.read_excel('../../../ESVT_data_allSubjects28Oct2022.xlsx')

all_ids = all_data['id_in_session'].unique()

df = pd.DataFrame(columns=['subject_id','ro', 'init_M', 'M', 'init_alpha', 'alpha', 'noise', 'likelihood'])


for id in all_ids:
    logging.warning(f"****************** Starting with id: {id} *************************")
    try:
        # if str(id) not in payload:
        #     continue

        data = all_data.loc[(all_data['id_in_session'] == id)  & (all_data['two_options'] == 1)]
        data.reset_index(drop=True, inplace=True)
        ro = data['alpha'].unique()[0]
        print(f"start for id={id}, ro={ro}")
        data_Pareto = data[data['pareto'] == 1]
        data_Uniform = data[data['pareto'] == 0]
        
        if len(data_Pareto) != len(data_Uniform) or len(data_Pareto) == 0 or len(data_Uniform) == 0:
            print(f"unequal lengths for id {id}")
            continue
        
        # Uniform
        data_Uniform.reset_index(drop=True, inplace=True)
        option1_high = np.array(data_Uniform['option1_high'])
        option1_low = np.array(data_Uniform['option1_low'])
        option2_high = np.array(data_Uniform['option2_high'])
        option2_low = np.array(data_Uniform['option2_low'])
        # ro = np.array(data_Uniform['alpha'])
        actual_chosen_options_data = np.array([0 if x==1 else 1 for x in list(data_Uniform['chosen_option'])])

        option1_high = np.power(option1_high, ro)
        option1_low = np.power(option1_low, ro)
        option2_high = np.power(option2_high, ro)
        option2_low = np.power(option2_low, ro)

        num_trials = len(data_Uniform)
        
        init_alpha_list = []
        init_M_list = []
        alpha_list = []
        noise_list = []
        M_list = []
        likelihood_list = []
        for i in range(10):
            M_bound = np.max(np.concatenate((option1_high, option1_low, option2_high, option2_low)))
            init_M = np.median(np.concatenate((option1_high, option1_low, option2_high, option2_low)))
            init_alpha = round(random.uniform(0.1, 5), 5)
            init_noise = 0.03
            initial_parameters = [init_alpha, init_M, init_noise]

            init_alpha_list.append(init_alpha)
            init_M_list.append(init_M)
            # Optimization using Nelder-Mead algorithm 
            result = minimize(obj_function, initial_parameters, (option1_high, option1_low, option2_high, option2_low, actual_chosen_options_data), method='Nelder-Mead', options={"disp": True, "maxiter": 1000}, tol=0.5, bounds=[(0.1, 5), (0, M_bound), (0, None)])

            optimized_parameters = result.x
            alpha_estimate = optimized_parameters[0]
            M_estimate = optimized_parameters[1]
            noise_estimate = optimized_parameters[2]
              
            likelihood_mean = obj_function([alpha_estimate, M_estimate, noise_estimate], option1_high, option1_low, option2_high, option2_low, actual_chosen_options_data)
            alpha_list.append(alpha_estimate)
            M_list.append(M_estimate)
            noise_list.append(noise_estimate)
            likelihood_list.append(likelihood_mean)

        logging.warning("------------------------------------------------------------------------------ RESULTS ------------------------------------------------------------------------------")
        best_idx = np.argmin(np.array(likelihood_list))
        list_row = [id, ro, init_M_list[best_idx], M_list[best_idx], init_alpha_list[best_idx], alpha_list[best_idx], noise_list[best_idx], likelihood_list[best_idx]]
        df.loc[len(df)] = list_row
        df.to_csv("Estimates_Uniform.csv", index = False)
        del data
    except Exception as e:
        logging.error(f"!!!!!!!!!!Error in id={id} : {str(e)}")

        
    logging.warning(f"****************** Done with id: {id} *************************")

