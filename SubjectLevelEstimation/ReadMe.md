The code for the estimations was run on NYU Langone's high performance computing cluster which utilizes GPU. The python scripts were executed using a batch script on Slurm. The code utilizes both Python libraries and CUDA for efficient computation with PyTorch.

We have used the python library scipy for the purpose of subject-level parameter estimation using maximum likelihood estimation (MLE). We specificlly minimize the neg-log-likelihood using scipy.optimize.minimize and utilize the Nelder-Mead Algorithm(Gao, F. and Han, L. Implementing the Nelder-Mead simplex algorithm with adaptive parameters. 2012. Computational Optimization and Applications. 51:1, pp. 259-277). For Robustness, we max-iteration limit of 1,000 and a stopping criteria of 0.5 tolerance. We initialized to the distributions’ medians, delta was initialized at 0.03, matching the sample-level pooled estimate. For the parameter, we took ten random initializations in the range {0.1,5} with a precision of 5. For calculating the likelihoods, in each of the 320 trials, we generated 10,000 samples with randomly drawn Gaussian noise. 

Requirements:

1. Python and Python Libraries
Please ensure you are using Python 3.8 or above. You can check your Python version with the following command: python --version
Please ensure the following libraries are installed in your environment:
pandas: For data manipulation and loading.
torch: PyTorch for GPU acceleration and tensor computations.
numpy: For numerical operations.
scipy: For the minimize function used in parameter optimization.
logging: To output logs and errors during execution.
They can be installed using the following pip command: pip install pandas torch numpy scipy

3. GPU Requirements
The program utilizes CUDA-enabled GPUs for faster computation with PyTorch. 
Please ensure that you have:
A compatible NVIDIA GPU with CUDA support.
CUDA drivers and PyTorch with CUDA installed.
Note: The code will work without GPUs, but will take very long for the estimations.
The computational requirements for this work were supported in part by the NYU Langone High Performance Computing (HPC) Core's resources.

The input data is expected to be in an Excel file named ESVT_data_allSubjects.xlsx, located in a relative path ../../

The program outputs CSV files containing the estimated parameters for each subject.

To run the program, use Python: python <script_name.py>
Note: We used Slurm jobs to submit the batch jobs on our HPC cluster. A sample job can be found in the file 'rscript.src'.
