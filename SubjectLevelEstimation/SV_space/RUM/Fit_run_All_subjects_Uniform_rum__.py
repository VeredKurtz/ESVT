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

# Objective Function
def obj_function_rum(alpha, option1_high, option1_low, option2_high, option2_low, actual_chosen_option):

    # Parameter Extraction
    delta = 0.68
    num_simulations_per_trial = 10000
    num_trials = len(option1_high)
    actual_chosen_option = torch.from_numpy(actual_chosen_option).to('cuda').float()
    option1_high = torch.from_numpy((np.repeat(option1_high, num_simulations_per_trial)).reshape(num_trials, num_simulations_per_trial)).to('cuda').float()
    option1_low = torch.from_numpy((np.repeat(option1_low, num_simulations_per_trial)).reshape(num_trials, num_simulations_per_trial)).to('cuda').float()
    option2_high = torch.from_numpy((np.repeat(option2_high, num_simulations_per_trial)).reshape(num_trials, num_simulations_per_trial)).to('cuda').float()
    option2_low = torch.from_numpy((np.repeat(option2_low, num_simulations_per_trial)).reshape(num_trials, num_simulations_per_trial)).to('cuda').float()
    alpha = alpha[0]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    N_option1_high = (option1_high**alpha)  
    N_option1_low = (option1_low**alpha) 
    N_option2_high = (option2_high**alpha) 
    N_option2_low = (option2_low**alpha)

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

all_data = pd.read_excel('../../ESVT_data_allSubjects28Oct2022.xlsx')

all_ids = all_data['id_in_session'].unique()

df = pd.DataFrame(columns=['subject_id','ro', 'mean_alpha', 'likelihood_mean'])


for id in all_ids:
    logging.warning(f"****************** Starting with id: {id} *************************")
    try:
        data = all_data.loc[(all_data['id_in_session'] == id)  & (all_data['two_options'] == 1)]
        data.reset_index(drop=True, inplace=True)
        ro = data['alpha'].unique()[0]
        print(f"start for id={id}, ro={ro}")
        data_Pareto = data[data['pareto'] == 1]
        data_Uniform = data[data['pareto'] == 0]
        
        if len(data_Pareto) != len(data_Uniform) or len(data_Pareto) == 0 or len(data_Uniform) == 0:
            print(f"unequal lengths for id {id}")
            continue
        
        del data

        # Uniform
        data_Uniform.reset_index(drop=True, inplace=True)
        option1_high = np.array(data_Uniform['option1_high'])
        option1_low = np.array(data_Uniform['option1_low'])
        option2_high = np.array(data_Uniform['option2_high'])
        option2_low = np.array(data_Uniform['option2_low'])
        actual_chosen_options_data = np.array([0 if x == 1 else 1 for x in list(data_Uniform['chosen_option'])])

        num_trials = len(data_Uniform)

        option1_high = np.power(option1_high, ro)
        option1_low = np.power(option1_low, ro)
        option2_high = np.power(option2_high, ro)
        option2_low = np.power(option2_low, ro)

        # Here population Estimates
        init_alpha = 0.5483
        result = minimize(obj_function_rum, [init_alpha], (option1_high, option1_low, option2_high, option2_low, actual_chosen_options_data), method='Nelder-Mead', options={"disp": True})

        optimized_parameters = result.x
        alpha_estimate = optimized_parameters[0]
        likelihood_mean = obj_function_rum([alpha_estimate], option1_high, option1_low, option2_high, option2_low, actual_chosen_options_data)
        logging.warning("------------------------------------------------------------------------------ RESULTS ------------------------------------------------------------------------------")
        list_row = [id, ro, alpha_estimate, likelihood_mean]
        df.loc[len(df)] = list_row
        df.to_csv("Estimates_Uniform_Rum.csv", index = False)

        
    except Exception as e:
        logging.error(f"!!!!!!!!!!Error in id={id} : {str(e)}")

    
    logging.warning(f"****************** Done with id: {id} *************************")

