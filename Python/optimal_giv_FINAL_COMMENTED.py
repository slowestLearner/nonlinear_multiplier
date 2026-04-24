####################################################################################
# THREADING AND ENVIRONMENT SETUP
####################################################################################
# Configures Numba and BLAS thread counts based on SLURM allocation or local CPUs.
# Must run before importing NumPy/SciPy so environment variables take effect.
####################################################################################

import os

def _int_env(name, default):
    try:
        return int(os.environ.get(name, default))
    except Exception:
        return default

def setup_threads():
    cpus = _int_env("SLURM_CPUS_PER_TASK", os.cpu_count() or 1)

    policy = os.environ.get("THREAD_POLICY", "numba").lower()
    if policy not in ("numba", "blas"):
        policy = "numba"

    if policy == "numba":
        numba_threads = cpus
        blas_threads  = 1
    else:
        numba_threads = 1
        blas_threads  = cpus

    os.environ["MKL_NUM_THREADS"]      = str(blas_threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(blas_threads)
    os.environ["OMP_NUM_THREADS"]      = str(blas_threads)
    os.environ["MKL_DYNAMIC"]          = "FALSE"
    os.environ["NUMBA_NUM_THREADS"]    = str(numba_threads)

    try:
        import numba as nb
        nb.set_num_threads(numba_threads)
    except Exception:
        pass

    try:
        import mkl
        mkl.set_num_threads(blas_threads)
    except Exception:
        pass

    try:
        from threadpoolctl import threadpool_limits
        globals()["_TP_CONTROLLER"] = threadpool_limits(limits=blas_threads,
                                                        user_api=["blas", "openmp"])
    except Exception:
        pass

    print(f"[threads] policy={policy}  numba={numba_threads}  blas={blas_threads}  "
          f"cpus={cpus}")

setup_threads()


####################################################################################
# IMPORTS
####################################################################################

import os
import glob
import pandas as pd
import pyreadr
import numpy as np
import time
import math
import sys
from scipy.optimize import fsolve, differential_evolution
import statsmodels.api as sm
import matplotlib.pyplot as plt
from copy import deepcopy
from functools import reduce, partial
import pickle
from pathlib import Path
from datetime import date, timedelta
import subprocess
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures  
import threading
from numba import jit, prange
from scipy.stats import qmc
from tqdm import tqdm
from scipy.optimize import least_squares, root
from hetero_pca import heteropca, predict, reconstruct, HeteroPCAModel
from hetero_pca_fast import heteropca_truly_fast
from scipy.optimize._numdiff import approx_derivative
from scipy import sparse
import numba as nb
from numba.typed import List as NbList

import warnings
warnings.filterwarnings("ignore")


####################################################################################
# WINSORIZATION FUNCTIONS
####################################################################################


def winsorize(df, series_name, lower_limit = .01, upper_limit = .99):
    """
    Winsorize a single column of a DataFrame in place.

    Inputs:
        df          : pd.DataFrame
        series_name : str, column name to winsorize
        lower_limit : float in [0,1], lower quantile cutoff (default 0.01)
        upper_limit : float in [0,1], upper quantile cutoff (default 0.99)

    Output:
        df : pd.DataFrame with df[series_name] clipped to [q_lower, q_upper]
    """
    lower_limit_value = df[series_name].quantile(lower_limit) 
    upper_limit_value = df[series_name].quantile(upper_limit) 

    df.loc[:, series_name] = df[series_name].clip(lower_limit_value, upper_limit_value)

    return df


def winsorize_np(a, limits=(0.01, 0.99), axis=None, in_place=False):
    """
    Winsorize array a by clipping values below/above the given quantile limits.
    
    Parameters
    ----------
    a : array-like
        Input data.
    limits : tuple of two floats in [0,1]
        (lower_quantile, upper_quantile). E.g. (0.01, 0.99) to clip the bottom/top 1%.
    axis : int or None
        Axis along which to compute quantiles. If None (default), flattens array.
    in_place : bool
        If True, clip in-place on `a` (must be a NumPy array).
        Otherwise returns a new array.
    
    Returns
    -------
    wins : ndarray
        Winsorized array.
    """
    a = np.asarray(a)
    # compute the lower and upper thresholds
    lower_q, upper_q = np.percentile(a, [limits[0]*100, limits[1]*100], axis=axis, keepdims=True)
    
    if in_place:
        # clip in-place
        a.clip(lower_q, upper_q, out=a)
        return a
    else:
        # return a new array
        return np.clip(a, lower_q, upper_q)
    
    


####################################################################################
# FILE READING FUNCTIONS
####################################################################################

file_cache = {}

try:
    num_workers = int(os.environ['SLURM_JOB_CPUS_PER_NODE']) - 1
except:
    num_workers = 4


def read_files_smart(file_paths):
    """
    Read .dta files with caching and parallel I/O.

    Inputs:
        file_paths : list of str, paths to .dta files

    Output:
        dict mapping file_path -> pd.DataFrame (with 'fdate_as_Q' column added)
    """
    cached = {f: file_cache[f].copy() for f in file_paths if f in file_cache}
    uncached = [f for f in file_paths if f not in file_cache]
    
    print(f"Files: {len(cached)} cached, {len(uncached)} reading...")
    
    if uncached:
        def read_and_cache(path):
            result = pd.read_stata(path)
            result['fdate_as_Q'] = result['fdate'].dt.to_period('Q')
            file_cache[path] = result.copy()
            return result
        
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            new_data = list(executor.map(read_and_cache, uncached))
        new_results = dict(zip(uncached, new_data))
    else:
        new_results = {}
    
    return {**cached, **new_results}


####################################################################################
# CUMULATIVE RETURNS
####################################################################################

def add_cum_return_calendar_fast(df, N):
    """
    Compute N-quarter cumulative return for each (permno, period).

    Requires N consecutive calendar quarters of data; gaps yield NaN.

    Inputs:
        df : pd.DataFrame with columns ['permno', 'fdate_as_Q', 'ret']
        N  : int, number of quarters to compound

    Output:
        df : pd.DataFrame with new column 'ret_cum_{N}q' (float)
    """
    df = df.copy()

    # Step 1: Parse quarters and sort
    df['period'] = pd.PeriodIndex(df['fdate_as_Q'], freq='Q')
    df = df.sort_values(['permno','period'])

    # Step 2: Build shifted returns ret_t, ret_{t-1}, ..., ret_{t-(N-1)}
    lag_rets = [ df.groupby('permno')['ret'].shift(i) for i in range(N) ]

    # Step 3: Get the period N-1 quarters back to verify calendar continuity
    df['period_lag'] = df.groupby('permno')['period'].shift(N-1)

    # Step 4: Compute ordinal quarter difference
    ord_p     = df['period'].astype(int)
    ord_plag  = df['period_lag'].astype(int)
    span      = ord_p - ord_plag

    # Step 5: Compound returns: prod_{k=0}^{N-1} (1 + ret_{t-k}) - 1
    prod = np.ones(len(df), dtype=float)
    for lag_r in lag_rets:
        prod *= (1 + lag_r)
    df[f'ret_cum_{N}q'] = prod - 1

    # Step 6: Set to NaN where the span is not exactly N-1 (i.e., missing quarters)
    df.loc[span != (N-1), f'ret_cum_{N}q'] = np.nan

    return df.drop(columns=['period_lag'])




####################################################################################
# RESIDUALIZATION: ORTHOGONALIZE ROWS W.R.T. OBSERVABLE/LATENT FACTORS
####################################################################################
# These functions project out stock characteristics eta_{n,t} from each row of a
# matrix (one row per investor), implementing the Frisch-Waugh-Lovell residualization
# described in Appendix C.2.1 of the paper.
####################################################################################

def _rowwise_residuals(Y, X):
    """
    OLS-residualize each row of Y against the common regressors X, handling NaNs.

    Inputs:
        Y : ndarray (I, N_t), one row per investor, one column per stock
        X : ndarray (N_t, K), regressors (stock characteristics / factors)

    Output:
        E : ndarray (I, N_t), residuals (NaN where Y was NaN)
    """
    N, T = Y.shape
    K = X.shape[1]

    E = np.full((N, T), np.nan)                         # (I, N_t)

    for i in prange(N):
        good = ~np.isnan(Y[i])
        Xi = X[good, :]                                  # (n_obs, K)
        yi = Y[i, good]                                  # (n_obs,)
        beta = np.linalg.pinv(Xi) @ yi                   # (K,)
        E[i, good] = yi - Xi @ beta                      # (n_obs,)

    return E


def residualize_rows(df: pd.DataFrame, F: np.ndarray, add_const=True):
    """
    Residualize every row of df against the factor matrix F.

    Inputs:
        df        : pd.DataFrame (I, N_t), rows = investors, columns = stocks
        F         : ndarray (K, N_t), factor loadings (each row is one factor)
        add_const : bool, whether to prepend a constant to the regressors

    Output:
        pd.DataFrame (I, N_t), residuals with same index/columns as df
    """
    X = F.T                                               # (N_t, K)
    if add_const:
        X = np.column_stack([np.ones(X.shape[0]), X])     # (N_t, K+1)

    E = _rowwise_residuals(df.values, X)                  # (I, N_t)
    return pd.DataFrame(E, index=df.index, columns=df.columns)



####################################################################################
# READ DISAGGREGATED MUTUAL FUND HOLDINGS DATA FROM .RDS FILES
####################################################################################


def read_rds_with_proper_dates(rds_file_path, year, local):
    """
    Read an RDS file via R subprocess, converting R date objects to Python dates.

    Inputs:
        rds_file_path : str, path to the .RDS file
        year          : int, used to create a unique temp CSV filename
        local         : bool, True if running locally (affects temp file paths)

    Output:
        pd.DataFrame with columns including 'rdate' and 'rdate_1' as Python dates,
        or None on failure.
    """
    
    # Create unique CSV filename with timestamp to avoid conflicts
    timestamp = int(time.time() * 1000)  # milliseconds


    if local:
        csv_filename = f'temp_dates_{year}_{timestamp}_{job_index}.csv' 

    else:
        csv_filename = f'/fs/scratch/PAS2771/temp_dates_{year}_{timestamp}_{job_index}.csv'
 

    # R script content
    r_script_content = f'''
    data <- readRDS('{rds_file_path}')
    
    # Convert dates to numeric (preserving original R date values)
    data$rdate_numeric <- as.numeric(data$rdate)
    data$rdate_1_numeric <- as.numeric(data$rdate_1)
    
    # Also keep original dates as character strings for verification
    data$rdate_char <- as.character(data$rdate)
    data$rdate_1_char <- as.character(data$rdate_1)
    
    # Remove original date columns to avoid confusion
    data$rdate <- NULL
    data$rdate_1 <- NULL
    
    # Save to CSV
    write.csv(data, '{csv_filename}', row.names=FALSE)
    cat("Saved to {csv_filename}\\n")
    '''
    
    # Write R script to temporary file
    if local:
        r_script_filename = f'temp_script_{year}_{timestamp}.R'

    else:
        r_script_filename = f'/fs/scratch/PAS2771/temp_script_{year}_{timestamp}.R'
 
    with open(r_script_filename, 'w') as f:
        f.write(r_script_content)
    
    try:
        # Run R script
        print(f"Processing {rds_file_path} via R...")
        result = subprocess.run(['Rscript', r_script_filename], 
                              capture_output=True, text=True)
        
        if result.returncode != 0:
            print("R script error:")
            print(result.stderr)
            return None
        
        # Check if CSV was created
        if not os.path.exists(csv_filename):
            print(f"CSV file {csv_filename} was not created")
            return None
        
        # Read CSV into pandas
        print(f"Reading {csv_filename} into Python...")
        df = pd.read_csv(csv_filename)
        
        # Convert numeric dates back to proper Python dates
        def numeric_to_date(numeric_val):
            if pd.isna(numeric_val):
                return None
            return date(1970, 1, 1) + timedelta(days=int(numeric_val))
        
        df['rdate'] = df['rdate_numeric'].apply(numeric_to_date)
        df['rdate_1'] = df['rdate_1_numeric'].apply(numeric_to_date)
        
        # Drop the temporary numeric and character columns
        df = df.drop(['rdate_numeric', 'rdate_1_numeric', 'rdate_char', 'rdate_1_char'], axis=1)
        
        print(f"Successfully loaded data with {len(df)} rows")
        print(f"Date range: {df['rdate'].min()} to {df['rdate'].max()}")
        
        return df
        
    except Exception as e:
        print(f"Error: {e}")
        return None
         
    finally:
        # Clean up temporary files
        try:
            if os.path.exists(r_script_filename):
                os.remove(r_script_filename)
                print(f"Deleted {r_script_filename}")
        except:
            pass 
             
        try:
            if os.path.exists(csv_filename):
                os.remove(csv_filename)
                print(f"Deleted {csv_filename}")
        except:
            pass


####################################################################################
# OUTPUT FILENAME CONSTRUCTION
####################################################################################






def get_output_filename(job_index, job_dict, local, outer_folder = None):
    """
    Build the pickle output path for the GMM estimation results.

    Inputs:
        job_index    : int, SLURM array task ID (indexes into job_dict)
        job_dict     : dict, maps job_index -> (job_file_index, num_factors)
        local        : bool, True if running locally
        outer_folder : unused (overridden internally)

    Output:
        str, path like 'HOLDINGS_ESTIMATION_RESULTS/numFactors_5/2010Q4.p'
    """
    job_file_index, num_factors = job_dict[job_index]

    outer_folder = 'HOLDINGS_ESTIMATION_RESULTS/'
    os.makedirs(outer_folder, exist_ok=True) 

    folder_name = '../../../Quarterly_Holdings_Full/' if local else 'Quarterly_Holdings_Full/'
    dta_files = sorted(glob.glob(os.path.join(folder_name, '*.dta')))
    target_file = dta_files[job_file_index + num_qtrs - 1]
    
    result = pd.read_stata(target_file)
    quarters = sorted(result['fdate'].dt.to_period('Q').unique())
    output_quarter = quarters[-1]

    os.makedirs(outer_folder, exist_ok=True) 

    subfolder = (f"numFactors_{num_factors}/")
    output_folder = outer_folder + subfolder 
    os.makedirs(output_folder, exist_ok=True)

    output_path = output_folder + str(output_quarter) + '.p'
     
    return output_path
 

def get_giv_output_filename(job_index, job_dict, local, outer_folder = None):
    """
    Build the .dta output path for the GIV (granular instrumental variable) results.

    Inputs:
        job_index    : int, SLURM array task ID
        job_dict     : dict, maps job_index -> (job_file_index, num_factors)
        local        : bool, True if running locally
        outer_folder : unused (overridden internally)

    Outputs:
        output_path     : str, e.g. '.../ESTIMATED_GIV_.../RedForm_GIV_quarter_2010Q4.dta'
        [output_quarter] : list of one Period, the last quarter in the estimation window
    """
    job_file_index, num_factors = job_dict[job_index]
    
    folder_name = '../../../Quarterly_Holdings_Full/' if local else 'Quarterly_Holdings_Full/'
    dta_files = sorted(glob.glob(os.path.join(folder_name, '*.dta'))) 
    
    target_file = dta_files[job_file_index + num_qtrs - 1]
    result = pd.read_stata(target_file)
    quarters = sorted(result['fdate'].dt.to_period('Q').unique())
    output_quarter = quarters[-1]

    scratch_dir = '/fs/scratch/PAS2771/'
    if local:
        scratch_dir = ''

    outer_folder = scratch_dir + ('ESTIMATED_GIV_numFactors_%s/' % (num_factors ) ) 
    os.makedirs(outer_folder, exist_ok=True) 

    subfolder = (f"numFactors_{num_factors}/")
    output_folder = outer_folder + subfolder 
    os.makedirs(output_folder, exist_ok=True)

    output_path = output_folder + 'RedForm_GIV_quarter_%s.dta' % output_quarter

    return output_path, [output_quarter]
 
  

 
############################################################################################################
# DETERMINE LOCAL vs. CLUSTER EXECUTION
############################################################################################################
DATA_DIR = Path('../../data/stocks/controls')
local = DATA_DIR.exists()
  
############################################################################################################
# READ STOCK CHARACTERISTICS (eta_{n,t} in the paper, i.e., observable controls)
# These serve as the observed component of eta_{n,t} in eq. (12) of the paper.
# Includes 13 return-predicting characteristics and Fama-French 12 industry dummies.
############################################################################################################
if local:
    DATA_DIR = Path('../../data/stocks/controls') 
    char_file = DATA_DIR / 'monthly_characteristics_not_lagged.RDS'
    ind_file  = DATA_DIR / 'ff12_industries_zero_mean.RDS'
else:
    char_file =  'monthly_characteristics_not_lagged.RDS' 
    ind_file  =   'ff12_industries_zero_mean.RDS'
   
# Step 1: Read firm characteristics
controls_df = pyreadr.read_r(char_file)[None]
vv_char = [c for c in controls_df.columns if c not in ('yyyymm', 'permno')]

# Step 2: Read industry classifications (FF12 dummies)
industries_df = pyreadr.read_r(ind_file)[None]
vv_ind = [c for c in industries_df.columns if c not in ('yyyymm', 'permno')]

controls_df = controls_df.merge(industries_df, on=['yyyymm', 'permno'], how = 'outer')
del industries_df

# Step 3: Convert monthly characteristics to quarterly, lagged by one quarter.
# Characteristics from month m are assigned to the NEXT quarter so they are
# available (pre-determined) when we observe holdings changes in that quarter.
controls_df = controls_df[['yyyymm', 'permno'] + vv_char + vv_ind]
controls_df = controls_df.dropna()
controls_df['yyyymm'] = controls_df['yyyymm'].astype(int)

dt = pd.to_datetime(controls_df['yyyymm'].astype(str) + '01', format='%Y%m%d')
q_next = dt.dt.to_period('Q') + 1 

controls_df['fdate_to_merge_on'] = q_next
controls_df = controls_df.sort_values( [ 'permno','fdate_to_merge_on' ] )
controls_df = controls_df.groupby( [ 'permno','fdate_to_merge_on' ] ).last().reset_index()

####################################################################################
# HOLDINGS DATA FOLDER
####################################################################################

if local: 
    folder_name = '../../../Quarterly_Holdings_Full/'
else:
    folder_name = 'Quarterly_Holdings_Full/'

dta_files = sorted(glob.glob(os.path.join(folder_name, '*.dta')))
data_frames = []

####################################################################################
# ESTIMATION SPECIFICATION PARAMETERS
####################################################################################
# Each SLURM array task runs one (job_file_index, num_factors) pair.
# job_file_index selects the starting quarter; num_factors selects the number
# of latent factors in eta_{n,t}:
#   -1 = no observed or latent factors
#    0 = only observed factors
#   K>0 = observed + K latent factors (estimated via deflated heteroskedastic PCA)
####################################################################################

job_dict = {}      
count = 0   
active_share_cutoff = -1  
max_ownership_share_threshold =  None

num_lagged_returns = 4              # L in the paper (number of lagged quarters)
num_qtrs = 4                        # rolling window length (in quarters)
wins_threshold = 0.0
prop_of_investors_needed = 20       # minimum number of investors per stock (Appendix C.2.2)
prop_of_stocks_needed = 60          # minimum number of stocks per investor (Appendix C.2.2)
char_list = ['active_share']        # investor characteristics Xi,t in eq. (18)
stock_chars_list = []

for num_factors in [ 5, -1, 0, 1, 3,  7, 10, 13, 15]:

    for job_file_index in range(3, len(dta_files) - 4):

        job_dict[count] = [ job_file_index , num_factors]
        count = count + 1

####################################################################################
# GET SPECIFICATION FOR THIS SLURM JOB
####################################################################################
if local:  
    job_index  = int(sys.argv[1])
else:
    job_index =   int(os.environ['SLURM_ARRAY_TASK_ID']) 

print('job_index: %s' % job_index) 

job_file_index , num_factors = job_dict[job_index] 

print('-----------------------')
print(job_file_index , num_factors)
print('-----------------------') 

####################################################################################
# BUILD OUTPUT PATHS
####################################################################################

output_path = get_output_filename(job_index, job_dict, local, None)
giv_output_path, output_quarters = get_giv_output_filename(job_index, job_dict, local, None)

print('OUTPUT QUARTERS: ', output_quarters)

outer_folder = output_path.split('/')[-2]


####################################################################################
# INITIALIZE ACCUMULATORS FOR ALL QUARTERS IN THE ROLLING WINDOW
####################################################################################
# Each list below will have length num_qtrs (one entry per quarter in the window).
# Each entry is itself a list of length (num_lagged_returns + 1), one per lag.
qmat_pd_list = []                   # Delta_q matrices (I x N_t DataFrames)
Cpts_pd_list = []                   # Delta_p * X_i matrices by investor type
Smat_pd_list = []                   # lagged ownership share matrices S_{i,n,t-1}
type_pd_list = []                   # investor type labels (e.g. bank, advisor, etc.)
Cpts_interaction_pd_list = []       # Delta_p * |P_tilde| * X_i interaction matrices
Cpts_lag_pd_list = []               # lagged cumulative return matrices
Cts_pd_list = []                    # X_i indicator matrices (used in moment conditions)

num_quaters = num_qtrs
have_data_for_non_active_mfs = False

####################################################################################
# READ QUARTERLY RETURNS AND COMPUTE CUMULATIVE RETURNS
####################################################################################

if local:
    returns_df = pyreadr.read_r( '../../data/stocks/prices/quarterly_return.RDS' )[None]

else:
    returns_df = pyreadr.read_r( 'quarterly_return.RDS' )[None]

returns_df['fdate_as_Q'] = pd.to_datetime(returns_df['yyyymm'].astype(int).astype(str), format = '%Y%m').dt.to_period('Q')
returns_df = returns_df.sort_values(['permno', 'fdate_as_Q'])

# Compute L-quarter cumulative returns (used to construct P_tilde_{n,t} in eq. (12))
returns_df = add_cum_return_calendar_fast( returns_df , N = num_lagged_returns )

####################################################################################
# READ INVESTOR CHARACTERISTICS (X_{i,t} in eq. (18): active share)
####################################################################################
# Active share is lagged by (L+1) quarters to avoid correlation with demand shocks
# from quarters t-L to t, as described in Section 5.1.3.
# We use the 4-quarter rolling mean, further shifted by (num_lagged_returns + 1).
####################################################################################

investor_chars_df = pd.read_csv('investor_ln_aum_inv_turnover_activeShare.csv')
investor_chars_df['fdate'] = pd.to_datetime(investor_chars_df['fdate'])
investor_chars_df = investor_chars_df.sort_values(['mgrno', 'fdate'])

all_fdates = investor_chars_df['fdate'].unique()
all_mgrnos = investor_chars_df['mgrno'].unique()
complete_grid = pd.MultiIndex.from_product([all_mgrnos, all_fdates], names=['mgrno', 'fdate']).to_frame(index=False)

investor_chars_df = pd.merge(complete_grid, investor_chars_df, on=['mgrno', 'fdate'], how='left')
investor_chars_df = investor_chars_df.sort_values(['mgrno', 'fdate'])
investor_chars_df = investor_chars_df[pd.notnull(investor_chars_df['fdate'])]

for col in char_list:
    if col in investor_chars_df.columns:
        investor_chars_df.loc[:, col] = investor_chars_df.groupby('mgrno')[col].transform(lambda x: x.rolling(window=4, min_periods=4).mean())
        investor_chars_df.loc[:, col] = investor_chars_df.groupby('mgrno')[col].shift( num_lagged_returns + 1) 

stock_chars_to_use_as_elasticity_interactions_df = controls_df[['permno', 'fdate_to_merge_on'] + stock_chars_list]

for col in stock_chars_list:
    if col in stock_chars_to_use_as_elasticity_interactions_df.columns:
        stock_chars_to_use_as_elasticity_interactions_df.loc[:, col] = stock_chars_to_use_as_elasticity_interactions_df.groupby('permno')[col].transform(lambda x: x.rolling(window=4, min_periods=4).mean())
        # Shift col forward by num_lagged_returns + 1 quarters within each mgrno
        stock_chars_to_use_as_elasticity_interactions_df.loc[:, col] = stock_chars_to_use_as_elasticity_interactions_df.groupby('permno')[col].shift( num_lagged_returns + 1) 


char_dict = {char : [] for char in char_list}


####################################################################################
# MAIN PROCESSING LOOP: READ AND PROCESS 13F HOLDINGS FOR EACH QUARTER
####################################################################################
# For each quarter in the rolling window, and for each lag l = 0, ..., L:
#   - Read 13F holdings data
#   - Merge with returns to get Delta_p_{n,t-l}
#   - Compute Delta_q_{i,n,t-l} (Davis-Haltiwanger percentage change)
#   - Compute lagged ownership shares S_{i,n,t-1-l}
#   - Build matrices for the GMM moment conditions (C.4)-(C.5)
####################################################################################

file_index_count = 0
dict_of_bin_mapping_dfs = {}

for file_index in tqdm(range(job_file_index, job_file_index + num_qtrs), desc="Processing quarters"):

    # Per-quarter accumulators: each list has length (num_lagged_returns + 1)
    this_q_qmat_pd_list = []             # Delta_q DataFrames (I x N_t)
    this_q_Cpts_pd_list = []             # Delta_p * X_i by type (list of lists)
    this_q_Smat_pd_list = []             # S_{i,n,t-1} DataFrames (I x N_t)
    this_q_type_list = []                # investor type DataFrames (I x N_t)
    this_q_Cpts_interaction_pd_list = [] # Delta_p * |P_tilde| * X_i (list of lists)
    this_q_Cpts_lag_pd_list = []         # lagged cumulative returns DataFrames
    this_q_Cts_pd_list = []              # X_i indicator DataFrames (list of lists)

    q_list = []

    start_time_this_quarter = time.time()

    # Step 1: Read and cache data files for this quarter (current + L lags)
    files_needed = [dta_files[file_index + lag] for lag in range(0, -num_lagged_returns - 1, -1)]
    file_data = read_files_smart(files_needed)

    file_index_count = file_index_count + 1

    # Memory management: clear cache if it exceeds threshold
    cache_size_gb = sum(df.memory_usage(deep=True).sum() for df in file_cache.values()) / (1024**3)

    if   local and cache_size_gb > 12:
        file_cache = {}
    elif not local and cache_size_gb > 20:
        file_cache = {}

    # Step 2: Process each lag l = 0, -1, ..., -L
    # lag=0 is the current quarter; lag=-l corresponds to l quarters back
    for lag in range(0, -num_lagged_returns -1, -1):

        file = dta_files[file_index + lag]
        result = file_data[file]
        result['fdate_as_Q'] = result['fdate'].dt.to_period('Q')

        quarters = sorted(result['fdate_as_Q'].unique())

        if max(quarters) < controls_df['fdate_to_merge_on'].min():
            sys.exit(0)

        if file_index == job_file_index + num_qtrs - 1 and lag == 0:
            output_quarters = quarters

        q_list.append(quarters[-1])

        if lag == 0:
            lag_0_quarter = quarters[-1]

        # Step 3: Merge with returns data (Delta_p_{n,t} and cumulative returns)
        result = pd.merge(left=result, 
                        right=returns_df[['fdate_as_Q', 'permno', 'ret', 'ret_cum_%sq' % num_lagged_returns]], 
                        how='left', on=['fdate_as_Q', 'permno'])

        # Step 4: Split into two quarters (Q1 = earlier, Q2 = later) for computing Delta_q
        result_q1 = result[result['fdate_as_Q'] == quarters[0]]
        result_q1 = result_q1.set_index(['mgrno', 'permno'])
        result_q1 = result_q1.rename(columns = {c : c + '_Q1' for c in result_q1.columns if c not in ['permno', 'mgrno']} )



        result_q2 = result[result['fdate_as_Q'] == quarters[1]]
        result_q2 = result_q2.set_index(['mgrno', 'permno'])
        result_q2 = result_q2.rename(columns = {c : c + '_Q2' for c in result_q2.columns if c not in ['permno', 'mgrno']} )

        # Merge the two quarters; missing positions are filled with 0 shares
        merged_result = pd.concat([result_q2 , result_q1], axis = 1)

        merged_result.loc[ pd.isnull(merged_result['shares_Q2']), 'shares_Q2'] = 0 
        merged_result.loc[ pd.isnull(merged_result['shares_Q1']), 'shares_Q1'] = 0

        # Flag extensive-margin entries (new positions or exits)
        merged_result['lag_position_0'] = 0
        merged_result.loc[ merged_result['shares_Q1']== 0 , 'lag_position_0'] = 1
        merged_result['current_position_0'] = 0
        merged_result.loc[ merged_result['shares_Q2']== 0 , 'current_position_0'] = 1

        merged_result.loc[ ( merged_result['shares_Q2'] > 0 ) & ( pd.isnull(merged_result['shares_Q1']) ) , 'shares_Q1'] = 0

        # Step 5: Merge investor-level characteristics (active share)

        merged_result = merged_result.reset_index()

        # All interaction characteristics use the lag-0 quarter (current quarter)
        sub_investor_chars_df = investor_chars_df[investor_chars_df['fdate'] == lag_0_quarter.to_timestamp(how = 'start')]
        sub_investor_chars_df = sub_investor_chars_df[['mgrno'] + char_list]

        for char in char_list:
            sub_investor_chars_df = winsorize(sub_investor_chars_df, char, 0.01, 0.99)

        merged_result = pd.merge(left = merged_result, right = sub_investor_chars_df, how = 'left', left_on = [ 'mgrno'], right_on = ['mgrno'], validate = 'm:1')


        for char in char_list:
            merged_result = merged_result[pd.notnull(merged_result[char])]


        sub_stock_chars_to_use_as_elasticity_interactions_df = stock_chars_to_use_as_elasticity_interactions_df.loc[stock_chars_to_use_as_elasticity_interactions_df['fdate_to_merge_on'] == lag_0_quarter] # Yes, all interaction characteristics should be from lag 0 quarter
        sub_stock_chars_to_use_as_elasticity_interactions_df = sub_stock_chars_to_use_as_elasticity_interactions_df.drop(columns=['fdate_to_merge_on'])

        merged_result = pd.merge(left=merged_result, 
                                right=sub_stock_chars_to_use_as_elasticity_interactions_df, 
                                how='left', 
                                on=['permno'],
                                validate='m:1')


        for char in stock_chars_list:
            merged_result = merged_result[pd.notnull(merged_result[char])]


        merged_result = merged_result.set_index(['mgrno', 'permno'])

        # Step 6: Apply intensive margin filter (restrict to observations where
        # investor holds positive quantities in both quarters, as in Section 5.1.3)
        merged_result = merged_result[merged_result['lag_position_0'] == 0]
        merged_result = merged_result[merged_result['current_position_0'] == 0]

        merged_result['type_to_use'] = 'all'

        # Step 7: Enforce consistent investor/stock samples across lags
        if lag == 0:

            mgr_list_from_base_quarter = merged_result.index.get_level_values('mgrno').unique()


        else:
            merged_result = merged_result.reset_index()


            merged_result = merged_result[merged_result['mgrno'].isin(mgr_list_from_base_quarter)]

            merged_result = merged_result[merged_result['permno'].isin(stock_list)]
            merged_result = merged_result[pd.notnull(merged_result['mgrno'])]
            merged_result = merged_result.set_index(['mgrno', 'permno'])

        # Aggregate holdings at the investor-stock level
        bin_delta_q = merged_result.groupby(['mgrno', 'permno'])[['shares_Q2', 'shares_Q1']].sum()
        bin_delta_q = bin_delta_q.reset_index()
        bin_delta_q.columns = ['mgrno', 'permno', 'bin_shares_Q2', 'bin_shares_Q1']

        merged_result = merged_result.reset_index()
        merged_result = pd.merge(left=merged_result, 
                                right=bin_delta_q.reset_index(), 
                                how='left', 
                                on=['mgrno', 'permno'],
                                validate='m:1')

        # Step 8: Compute and demean investor-level active share
        # Active share is demeaned cross-sectionally within each quarter,
        # so X_{i,t} = [1, ActiveShare_i - mean(ActiveShare)] as in eq. (18).
        if 'active_share' in char_list:
            bin_avg_active_share = merged_result.groupby('mgrno')[['active_share']].mean().reset_index()
            bin_avg_active_share.columns = ['mgrno', 'active_share']

            bin_avg_active_share.loc[:, 'active_share'] =  ( bin_avg_active_share.loc[:, 'active_share'] - bin_avg_active_share['active_share'].mean())

            merged_result = merged_result.drop(columns=['active_share'])
            merged_result = pd.merge(left=merged_result, 
                                    right=bin_avg_active_share,
                                    how='left', 
                                    on='mgrno',
                                    validate='m:1')



        # Step 9: Calculate Delta_q_{i,n,t} using Davis-Haltiwanger (1992) formula
        # Delta_q = (Q_hat_t - Q_hat_{t-1}) / (0.5 * (Q_hat_t + Q_hat_{t-1}))
        # where Q_hat_{t-1} = dollar holdings, Q_hat_t = dollar holdings / (1 + R^X)
        # This maps to [-2, 2], effectively winsorizing large percentage changes.
        # See footnote 15 of the paper.


        merged_result['within_bin_ownership_share_Q1'] = merged_result.loc[:, 'shares_Q1'] / merged_result.loc[:, 'bin_shares_Q1']

        # Ex-dividend return R^X_{n,t}
        merged_result['retx'] = ( np.exp(merged_result['LNprc_Q2']) * np.exp(merged_result['LNshrout_Q2']) - np.exp(merged_result['LNprc_Q1']) * np.exp(merged_result['LNshrout_Q1']) ) / ( np.exp(merged_result['LNprc_Q1']) * np.exp(merged_result['LNshrout_Q1']))

        # Q_hat_{t-1} = H_{i,n,t-1} (lagged dollar holdings)
        lagged_dollar_holdings = merged_result['shares_Q1'] *  np.exp(merged_result['LNprc_Q1']) 
        # Q_hat_t = H_{i,n,t} / (1 + R^X) (current dollar holdings deflated by ex-div return)
        contemp_dollar_holdings_deflated = merged_result['shares_Q2'] *  np.exp(merged_result['LNprc_Q2'])  / ( 1 + merged_result['retx']  )

        # Davis-Haltiwanger percentage change: Delta_q_{i,n,t}
        delta_q = ( contemp_dollar_holdings_deflated - lagged_dollar_holdings ) / ( .5 * ( contemp_dollar_holdings_deflated + lagged_dollar_holdings ) )
    
        merged_result['delta_q_inv_level'] = delta_q
        merged_result = winsorize(merged_result, 'delta_q_inv_level', wins_threshold, 1 - wins_threshold)

        # Remove investor-quarter fixed effects from Delta_q (Section 5.1.3)
        merged_result['delta_q_inv_level'] = merged_result['delta_q_inv_level'] - \
                                        merged_result.groupby(['mgrno'])['delta_q_inv_level'].transform('mean')

        merged_result['delta_q'] = merged_result['delta_q_inv_level']

        merged_result.loc[:, 'shares_Q1'] = merged_result.loc[:, 'bin_shares_Q1']
        merged_result.loc[:, 'shares_Q2'] = merged_result.loc[:, 'bin_shares_Q2']
 
        merged_result = merged_result.groupby(['mgrno', 'permno']).first()

        # S_{i,n,t-1} = lagged ownership share (proportion of shares outstanding)
        merged_result['lagged_ownership_share'] = merged_result['shares_Q1'] / np.exp(merged_result['LNshrout_Q1'])

        # Step 10: Select key variables for matrix construction

        merged_result = merged_result[['fdate_Q2',  'delta_q',  'type_Q2', 'aum_Q2', 'current_position_0', 
                                        'lag_position_0', 'retx', 'ret_Q2', 'ret_Q1', 'lagged_ownership_share', 
                                        'LNme_Q1', 'LNbe_Q1', 'ret_cum_%sq_Q1' % num_lagged_returns , 'ret_cum_%sq_Q2' % num_lagged_returns,
                                        'type_to_use'] + char_list + stock_chars_list] 

    
        merged_result.loc[:, 'delta_q'] = merged_result.loc[:, 'delta_q'] * 100
        merged_result.loc[:, 'ret_to_use'] = merged_result.loc[:, 'ret_Q2'] * 100
        merged_result.loc[:, 'ret_to_use_lag'] = merged_result.loc[:, 'ret_cum_%sq_Q1' % num_lagged_returns] * 100


        merged_result = merged_result[pd.notnull(merged_result['ret_to_use'])]
        merged_result = merged_result[pd.notnull(merged_result['delta_q'])]
        merged_result = merged_result[pd.notnull(merged_result['lagged_ownership_share'])]
        merged_result = merged_result[pd.notnull(merged_result['ret_to_use_lag'])]


        merged_result = merged_result.rename(columns = {'fdate_Q2' : 'fdate', 'type_Q2' : 'type'} )
        merged_result = merged_result.reset_index()


        # Apply minimum holdings and investor filters (Appendix C.2.2):
        # keep investors with >= prop_of_stocks_needed stocks
        # keep stocks held by >= prop_of_investors_needed investors
        merged_result.loc[:, 'num_holdings_per_investor'] = merged_result.groupby('mgrno')['delta_q'].transform(lambda x: x.count())
        merged_result = merged_result[merged_result['num_holdings_per_investor'] >= prop_of_stocks_needed]

        merged_result.loc[:, 'num_investors_per_stock'] = merged_result.groupby('permno')['delta_q'].transform(lambda x: x.count())
        merged_result = merged_result[merged_result['num_investors_per_stock'] >=    prop_of_investors_needed]

        merged_result = winsorize(merged_result, 'ret_to_use', .01, .99)

        mgrno_list = list(merged_result['mgrno'].unique())
        permno_list = list(merged_result['permno'].unique())

        if lag == 0:
            stock_list = sorted(permno_list)

        # Step 11: Pivot from long format to (I x N_t) matrices for GMM
        q_avg_inv_q_fe = merged_result.loc[:, 'delta_q'].mean()

        Smat_pd = merged_result.reset_index().pivot(index='mgrno', columns='permno', values='lagged_ownership_share')  # (I, N_t)
        this_q_Smat_pd_list.append(Smat_pd)

        qmat_pd = merged_result.reset_index().pivot(index='mgrno', columns='permno', values='delta_q')                # (I, N_t)
        this_q_qmat_pd_list.append(qmat_pd)

        type_pd = merged_result.reset_index().pivot(index='mgrno', columns='permno', values='type')                   # (I, N_t)
        this_q_type_list.append(type_pd)

        type_to_use_pd = merged_result.reset_index().pivot(index='mgrno', columns='permno', values='type_to_use')     # (I, N_t)

        # Delta_p_{n,t-l} matrix (contemporaneous return for this lag)
        Cpts_pd = merged_result.reset_index().pivot(index='mgrno', columns='permno', values='ret_to_use')             # (I, N_t)

        # Cumulative lagged return matrix (used to construct |P_tilde_{n,t}|)
        Cpts_lag_pd = merged_result.reset_index().pivot(index='mgrno', columns='permno', values='ret_to_use_lag')      # (I, N_t)
        this_q_Cpts_lag_pd_list.append(Cpts_lag_pd)

        # |P_tilde_{n,t}| = |sum_{l=1}^L Delta_p_{n,t-l}| demeaned across stocks
        matrix_for_Cts_pd = np.abs(Cpts_lag_pd)                                                                        # (I, N_t)

        matrix_for_Cts_pd = matrix_for_Cts_pd.sub( matrix_for_Cts_pd.mean(1).values, axis = 0)  # demean across stocks for each investor

        # Interaction term: Delta_p_{n,t} * P_tilde_{n,t} (used for zeta_2 identification)
        Cpts_interaction_pd = Cpts_pd * matrix_for_Cts_pd.values
        d =  pd.melt(Cpts_interaction_pd, ignore_index = False).reset_index().set_index(['mgrno', 'permno'])
        d = winsorize(d, 'value', 0.01, 0.99)
        Cpts_interaction_pd = d.reset_index().pivot(index='mgrno', columns='permno', values='value')



        # Step 12: Build per-type and per-characteristic moment condition matrices
        # For each investor type t and characteristic k, create:
        #   Cpts * I(type=t) * X_{i,k}  (used for zeta_1 moments, eq. C.4)
        #   Cpts_interaction * I(type=t) * X_{i,k}  (used for zeta_2 moments, eq. C.5)
        #   I(type=t) * X_{i,k}  (the C_{i,t} weighting matrices)
        investor_type_list = sorted(list(merged_result['type'].unique()))
        investor_type_to_use_list = sorted(list(merged_result['type_to_use'].unique()))

        Cpts_times_type_pd_list = [ Cpts_pd * (type_to_use_pd == t).values  for t in investor_type_to_use_list]
        Cpts_interaction_times_type_pd_list = [ Cpts_interaction_pd * (type_to_use_pd == t).values  for t in investor_type_to_use_list]


        ones_df = deepcopy(Cpts_pd)
        ones_df.loc[:, :] = 1
        
        all_Cts_pd_list = [ ones_df *  (type_to_use_pd == t).values   for t in investor_type_to_use_list]  

        # Extend lists with characteristic interactions (e.g., active_share * type indicator)
        for char in char_list + stock_chars_list:
            char_pd = merged_result.reset_index().pivot(index='mgrno', columns='permno', values=char)
            char_pd = char_pd[Cpts_pd.columns]


            Cpts_times_char_pd_list = [ Cpts_times_type_pd_list[i] * char_pd.values for i in range(len(investor_type_to_use_list))]
            Cpts_times_type_pd_list.extend(Cpts_times_char_pd_list)


            if char in char_list:
                Cpts_interaction_times_char_pd_list = [ Cpts_interaction_times_type_pd_list[i] * char_pd.values for i in range(len(investor_type_to_use_list))]
                Cpts_interaction_times_type_pd_list.extend(Cpts_interaction_times_char_pd_list)


            new_Cts_pd_list = [ all_Cts_pd_list[i] * char_pd.values for i in range(len(investor_type_to_use_list))]
            all_Cts_pd_list.extend(new_Cts_pd_list)


        this_q_Cpts_pd_list.append(Cpts_times_type_pd_list)
        this_q_Cpts_interaction_pd_list.append(Cpts_interaction_times_type_pd_list)
        this_q_Cts_pd_list.append(all_Cts_pd_list)

    # Step 13: Ensure consistent (I x N_t) dimensions across all lags
    # Find the intersection of investors and stocks present in every lag
    common_mgrno_list = sorted(list(this_q_qmat_pd_list[0].index))
    common_permno_list = sorted(list(this_q_qmat_pd_list[0].columns))


    for index in range(len(this_q_qmat_pd_list)):
        l = sorted(list(this_q_qmat_pd_list[index].index))
        common_mgrno_list = set(common_mgrno_list).intersection(l)
        common_mgrno_list = sorted(list(common_mgrno_list))

        l = sorted(list(this_q_qmat_pd_list[index].columns))
        this_quarter = q_list[index]
        sub_controls_df = controls_df[controls_df['fdate_to_merge_on'] == this_quarter]
        permnos_we_have_chars_for = list(sub_controls_df['permno'].unique())

        l = set(l).intersection(permnos_we_have_chars_for)
        l = sorted(list(l))

        common_permno_list = set(common_permno_list).intersection(l)
        common_permno_list = sorted(list(common_permno_list))

    # Subset all matrices to common (I x N_t) dimensions
    this_q_Smat_pd_list = [d.loc[common_mgrno_list, common_permno_list] for d in this_q_Smat_pd_list]
    this_q_type_list = [d.loc[common_mgrno_list, common_permno_list] for d in this_q_type_list]



    # Step 14: Residualize Delta_q and Delta_p w.r.t. stock characteristics eta_{n,t}
    # This implements the Frisch-Waugh-Lovell residualization from Appendix C.2.1:
    #   check_q_{i,n,t} = Delta_q - projection onto [observed chars, latent factors]
    #   check_p_{n,t}   = Delta_p - projection onto [observed chars, latent factors]
    # Doing so is equivalent to controlling for eta_{n,t} in the demand curve.
    for l in range(num_lagged_returns + 1):
        qmat_pd = this_q_qmat_pd_list[l]
        Cpts_times_type_pd_list = this_q_Cpts_pd_list[l]
        Cpts_interaction_times_type_pd_list = this_q_Cpts_interaction_pd_list[l]
        all_Cts_pd_list = this_q_Cts_pd_list[l] 
        
        # Step 14a: Remove investor-quarter and global fixed effects from Delta_q
        qmat_pd = qmat_pd.loc[common_mgrno_list, common_permno_list].sub(qmat_pd.loc[common_mgrno_list, common_permno_list].mean(1).values, axis = 0)
        qmat_pd = qmat_pd - np.nanmean(qmat_pd)

        qmat_pd = qmat_pd - np.nanmean(qmat_pd, axis = 1)[:, None]

        # Step 14b: Remove observable characteristics from Delta_q
        if num_factors > -1:
            this_quarter = q_list[l]
            sub_controls_df = controls_df[controls_df['fdate_to_merge_on'] == this_quarter]
            sub_controls_df = sub_controls_df[ sub_controls_df['permno'].isin(common_permno_list) ]

            V_obs_chars = sm.add_constant( sub_controls_df.iloc[:, 3:]).values.T  # (K_obs+1, N_t)

            qmat_pd = residualize_rows(qmat_pd, V_obs_chars, add_const=False)
            qmat_pd = qmat_pd - np.nanmean(qmat_pd)

        else:
            V_obs_chars = np.ones( ( 1 , len(common_permno_list) ) )    # (1, N_t) just a constant

        # Step 14c: Remove latent factors (deflated heteroskedastic PCA, Zhou & Chen 2025)
        # Estimates num_factors latent stock characteristics from the Delta_q panel
        qmat = qmat_pd.values                                          # (I, N_t)

        Smat_pd = this_q_Smat_pd_list[ l ]
        Smat = Smat_pd.values                                          # (I, N_t)

        exclude_mask = np.isnan(qmat) + np.isnan(Smat)
        Smat = Smat.copy()
        Smat[np.isnan(Smat)] = 0.0        
        mask     = ~exclude_mask                                        # (I, N_t) bool
        mf       = mask.astype(qmat.dtype)                              # (I, N_t) float

        if num_factors > 0:

            
            delta_q_for_pca = pd.melt(qmat_pd, ignore_index = False, value_name = 'delta_q').reset_index()

            wins_threshold_for_pca = 0.1
            delta_q_for_pca['delta_q_moreWins'] = winsorize(deepcopy(delta_q_for_pca), 'delta_q', wins_threshold_for_pca, 1 - wins_threshold_for_pca)['delta_q']
            delta_q_for_pca = delta_q_for_pca[delta_q_for_pca['permno'].isin( qmat_pd.columns)]
            

            X = pd.pivot_table(delta_q_for_pca, index = ['mgrno'], columns = 'permno', values = 'delta_q_moreWins')
            X = X.values

            model = heteropca_truly_fast(X, rank = num_factors, impute_method = 'zero', maxiter = 50000, verbose = False ) 

            V = predict(model, X, lambda_reg = 1e-8)              # (num_factors, N_t)
            qmat_with_factors_removed = residualize_rows(qmat_pd, V, add_const=False)  

            # Stack observed and latent factors: V = [V_obs; V_latent_demeaned]
            V = np.vstack( [ V_obs_chars , V - V.mean(axis=1, keepdims=True)] )  # (K_obs+1+num_factors, N_t)


        else:
            qmat_with_factors_removed = qmat
            V = V_obs_chars

        qmat_with_factors_removed =  qmat_with_factors_removed - np.nanmean(qmat_with_factors_removed, axis = 1)[:, None]

        qmat_pd = pd.DataFrame(qmat_with_factors_removed, index = qmat_pd.index, columns = qmat_pd.columns)  # (I, N_t)

        this_q_qmat_pd_list[l] = qmat_pd

        # Step 14d: Apply same residualization to Delta_p and Delta_p * P_tilde matrices

        char_mean_list = []
        for type_index in range(len(Cpts_times_type_pd_list)):
            Cpts_pd = Cpts_times_type_pd_list[type_index]
            Cts_pd = all_Cts_pd_list[type_index]


            Cpts_pd = Cpts_pd.loc[common_mgrno_list, common_permno_list]                    # (I, N_t)
            Cts_pd = Cts_pd.loc[common_mgrno_list, common_permno_list]                     # (I, N_t)

            if type_index > 0 and type_index < 1 + len(char_list):
                char_mean = deepcopy(Cts_pd.max(1).mean())
                char_mean_list.append(char_mean)

                Cts_pd = Cts_pd - char_mean
                Cpts_pd = Cpts_pd - char_mean * Cpts_times_type_pd_list[0]                

            Cpts_pd = Cpts_pd - np.nanmean(Cpts_pd)

            if num_factors > -1:
                Cpts_pd = residualize_rows(Cpts_pd, V, add_const = False)

            Cpts_pd = Cpts_pd - np.nanmean(Cpts_pd)


            Cpts_times_type_pd_list[type_index] = Cpts_pd
            all_Cts_pd_list[type_index] = Cts_pd


        for type_index in range(len(Cpts_interaction_times_type_pd_list)):
            Cpts_interaction_pd = Cpts_interaction_times_type_pd_list[type_index]

            Cpts_interaction_pd = Cpts_interaction_pd.loc[common_mgrno_list, common_permno_list] 

            if type_index > 0 and type_index < 1 + len(char_list):
                Cpts_interaction_pd = Cpts_interaction_pd - char_mean_list[type_index - 1] * Cpts_interaction_times_type_pd_list[0]


            Cpts_interaction_pd = Cpts_interaction_pd - np.nanmean(Cpts_interaction_pd)

            if num_factors > -1:
                Cpts_interaction_pd = residualize_rows(Cpts_interaction_pd, V, add_const = False)

            Cpts_interaction_pd = Cpts_interaction_pd - np.nanmean(Cpts_interaction_pd)


            Cpts_interaction_times_type_pd_list[type_index] = Cpts_interaction_pd


        this_q_Cpts_pd_list[l] = Cpts_times_type_pd_list
        this_q_Cpts_interaction_pd_list[l] = Cpts_interaction_times_type_pd_list
        this_q_Cts_pd_list[l] = all_Cts_pd_list

    # Step 15: Store this quarter's processed matrices into the master lists
    qmat_pd_list.append(this_q_qmat_pd_list)
    Cpts_pd_list.append(this_q_Cpts_pd_list)
    Smat_pd_list.append(this_q_Smat_pd_list)
    type_pd_list.append(this_q_type_list)
    Cpts_interaction_pd_list.append(this_q_Cpts_interaction_pd_list)
    Cts_pd_list.append(this_q_Cts_pd_list)



del file_cache


####################################################################################
# TENSOR CREATION: CONVERT QUARTERLY DATAFRAMES INTO 3D/4D NUMPY ARRAYS
####################################################################################
# Concatenate quarters within each lag to form tensors with dimensions:
#   (L+1, I, N)  for Delta_q, S, type
#   (L+1, I, N, P) for Delta_p-related moment condition matrices
# where L = num_lagged_returns, I = investors, N = total stocks across all
# quarters in the rolling window, P = number of moment condition types.
####################################################################################

def get_3d_tensor(df_list):
    """
    Stack quarterly DataFrames into a 3D array (one slice per lag).

    Inputs:
        df_list : list of length num_qtrs, each element is a list of length (L+1)
                  of DataFrames (I x N_t)

    Output:
        ndarray (L+1, I, N) where N = sum of N_t across quarters
    """
    mat_list = []
    for lag in range(num_lagged_returns + 1):
        l = []
        for qtr in range(num_quaters):
            l.append(df_list[qtr][lag])
        mat = pd.concat(l, axis = 1).values
        mat_list.append( mat )
    return np.array(mat_list)


    

def get_3d_tensor_Cpts(df_list, include_stock_chars = False):
    """
    Stack quarterly per-type DataFrames into a 4D array.

    Inputs:
        df_list            : nested list [quarter][lag][type_index] of DataFrames (I x N_t)
        include_stock_chars : bool, whether to include stock characteristic interactions

    Output:
        ndarray (L+1, I, N, P) where P = num_types * (1 + len(char_list) + ...) is the
        number of moment condition types per lag
    """
    mat_list = []
    for lag in range(num_lagged_returns + 1):
        this_lag_mat_list = []
        
        for type_index in range( len(investor_type_to_use_list) * (len(char_list) + int(include_stock_chars) * len(stock_chars_list) + 1) ):

            l = []
            for qtr in range(num_quaters):
                l.append(df_list[qtr][lag][type_index])
            
            mat = pd.concat(l, axis = 1).values[:, :, None]    # (I, N, 1)
            this_lag_mat_list.append( mat )
        
        mat = np.concatenate( this_lag_mat_list, axis = 2 )     # (I, N, P)
        mat_list.append(mat)
    
    return np.array(mat_list)                                    # (L+1, I, N, P)


    


####################################################################################
# CREATE MAIN DATA TENSORS FOR GMM ESTIMATION
####################################################################################

qmat = get_3d_tensor(qmat_pd_list)                                                     # (L+1, I, N) residualized Delta_q
Smat = get_3d_tensor(Smat_pd_list)                                                     # (L+1, I, N) lagged ownership shares S_{i,n,t-1}
type_tensor = get_3d_tensor(type_pd_list)                                               # (L+1, I, N) investor type labels

investor_type_list = sorted(list(merged_result['type'].unique()))

# Delta_p and Delta_p * P_tilde tensors, split by investor type x characteristic
Cpts = get_3d_tensor_Cpts(Cpts_pd_list, include_stock_chars = True)                     # (L+1, I, N, P1)
Cpts_interaction = get_3d_tensor_Cpts(Cpts_interaction_pd_list, include_stock_chars = False)  # (L+1, I, N, P2)

# Stack zeta_1 and zeta_2 moment condition matrices along the parameter axis
# First P1 columns identify zeta_1 params; next P2 identify zeta_2 params
Cpts =  np.concatenate( [ Cpts, 1 * Cpts_interaction ], axis = 3)                      # (L+1, I, N, P1+P2)

# Corresponding X_{i,t} indicator matrices (for weighting in moment conditions)
Cts_raw_term = get_3d_tensor_Cpts(Cts_pd_list, include_stock_chars = True)              # (L+1, I, N, P1)
Cts_raw_interaction_term = get_3d_tensor_Cpts(Cts_pd_list, include_stock_chars = False)  # (L+1, I, N, P2)
Cts =  np.concatenate( [ Cts_raw_term, Cts_raw_interaction_term ], axis = 3)            # (L+1, I, N, P1+P2)
Cts[np.isnan(Cts)] = 0.0

# Nmom_list[0] = number of zeta_1 moment conditions (for lag 0)
# Nmom_list[1] = number of zeta_2 moment conditions (for lags 1..L)
Nmom_list = [len(investor_type_to_use_list) * (len(char_list) + len(stock_chars_list) + 1) ,len(investor_type_to_use_list) * (len(char_list) + 1) ]
print(Nmom_list)

####################################################################################
# MISSING DATA HANDLING AND LEAVE-ONE-TYPE-OUT MATRIX
####################################################################################

exclude_mask = np.isnan(qmat) + np.isnan(Smat)                                         # (L+1, I, N) bool
Smat[np.isnan(Smat)] = 0.0

mask     = ~exclude_mask                                                                # (L+1, I, N) bool
mf       = mask.astype(qmat.dtype)                                                      # (L+1, I, N) float

# Leave-one-type-out matrix (Appendix C.2.2):
# neq_mat[i,j] = 1 if investors i and j are of DIFFERENT 13F institution types, 0 if same.
# GIV instruments z_{i,n,t} sum over j != i of a different type (eq. 15).
type_pd = pd.DataFrame(type_tensor[0])
type_vec = type_pd.T.ffill().T.iloc[:, -1]
vals = type_vec.to_numpy()

type_neq_mat = (vals[:, None] != vals).astype(int)                                      # (I, I)
neq_mat =  type_neq_mat 

unique_vals = np.unique(vals)
block_id = np.searchsorted(unique_vals, vals)                                           # (I,) integer block IDs

# GMM weighting matrix: inversely proportional to number of common observations
C = 1/mf[0].dot(mf[0].T) * neq_mat                                                    # (I, I)
C[ mf[0].dot(mf[0].T) == 0 ] = 0.0
Cdiag_col = np.diag(C)[:, None]                                                        # (I, 1)





####################################################################################
# GMM MOMENT CONDITION COMPUTATION
####################################################################################
# These functions compute the moment conditions (C.6)-(C.7) from Appendix C.2.
# The key objects are:
#   u_{i,n,t}(zeta) = check_q_{i,n,t} + zeta' * check_p_{n,t}   (residual demand shock, eq. C.3)
#   z_{i,n,t}(zeta) = sum_{j != i} S_{j,n,t-1} * u_{j,n,t}      (GIV instrument, eq. 15)
#   Z_tilde_{i,n,t} = |sum_l z_{i,n,t-l}| - E^CX[|sum_l z_{i,n,t-l}|]
#
# The moment conditions are:
#   (C.6): sum_i w_i * E^CX[ X_{i,t,k} * u_{i,n,t} * z_{i,n,t} ] = 0
#   (C.7): sum_i w_i * E^CX[ X_{i,t,k} * u_{i,n,t} * Z_tilde * z_{i,n,t} ] = 0
#
# where w_i are precision weights (inverse residual variance, eq. C.8).
####################################################################################
 
@nb.njit(parallel=True, fastmath=True)
def block_sum_numba_t(u, s, block_id, G):
    """
    Compute block-wise sums of u*s grouped by investor type, parallelized over stocks.

    Inputs:
        u        : ndarray (I, N), residual demand shocks
        s        : ndarray (I, N), ownership shares
        block_id : ndarray (I,), integer type label for each investor
        G        : int, number of distinct investor types

    Output:
        out : ndarray (G, N), sum of u*s within each type block
    """
    N, T = u.shape
    out = np.zeros((G, T), dtype=u.dtype)
    for t in nb.prange(T):
        for i in range(N):
            g = block_id[i]
            val = u[i, t] * s[i, t]
            if val != 0.0:
                out[g, t] += val
    return out

def build_CB_cache_numba(umat, Sm, block_id, use_diff=True):
    """
    Build the leave-one-type-out GIV z_{i,n,t} for each lag using block sums.

    z_{i,n,t} = sum_{j: type(j) != type(i)} S_{j,n,t-1} * u_{j,n,t}
              = (total across all types) - (sum within investor i's own type)

    Inputs:
        umat     : ndarray (L+1, I, N), residual demand shocks
        Sm       : ndarray (L+1, I, N), ownership shares (cleaned)
        block_id : ndarray (I,), integer type label for each investor
        use_diff : bool, if True compute leave-one-type-out (total - own block)

    Output:
        list of length L+1, each element is ndarray (I, N) = z_{i,n,t-l}
    """
    L, N, T = umat.shape
    G = int(block_id.max()) + 1
    cb = []
    for lag in range(L):
        blk_sum = block_sum_numba_t(umat[lag], Sm[lag], block_id, G)   # (G, N)
        if use_diff:
            total = blk_sum.sum(0)                                      # (N,)
            CB = total - blk_sum[block_id]                              # (I, N)
        else:
            CB = blk_sum[block_id] - umat[lag]*Sm[lag]
        cb.append(CB)
    return cb
    

@nb.njit(parallel=True, fastmath=True)
def build_umat_sm64(Cpts, qmat, Smat, zeta, mf):
    """
    Construct the residual demand shock matrix u_{i,n,t}(zeta) and cleaned S matrix.

    u_{i,n,t} = check_q_{i,n,t} + zeta' * check_p_{n,t}   (eq. C.3)

    Inputs:
        Cpts : ndarray (L+1, I, N, P), residualized Delta_p stacked by moment type
        qmat : ndarray (L+1, I, N), residualized Delta_q
        Smat : ndarray (L+1, I, N), lagged ownership shares
        zeta : ndarray (P,), current parameter guess
        mf   : ndarray (L+1, I, N), data availability mask (1.0 or 0.0)

    Outputs:
        umat : ndarray (L+1, I, N), residual demand shocks (0 where masked)
        Sm   : ndarray (L+1, I, N), cleaned ownership shares (0 where masked)
    """
    zeta64 = zeta.astype(np.float64)
    L, N, T, P = Cpts.shape
    umat = np.empty((L, N, T), dtype=np.float64)                       # (L+1, I, N)
    Sm   = np.empty_like(umat)                                          # (L+1, I, N)

    for l in nb.prange(L):
        for i in range(N):
            for t in range(T):
                if not mf[l, i, t]:
                    umat[l, i, t] = 0.0
                    Sm[l, i, t]   = 0.0
                    continue
                acc = qmat[l, i, t]
                for k in range(P):
                    acc += zeta64[k] * Cpts[l, i, t, k]     # now pure f64 math
                umat[l, i, t] = 0.0 if np.isnan(acc) else acc
                s = Smat[l, i, t]
                Sm[l, i, t]   = 0.0 if np.isnan(s) else s
    return umat, Sm


@nb.njit(parallel=True, fastmath=True)
def step7_numba(base0, base1, CB_all, Cts, mf, Nmom0, Nmom1, err_out):
    """
    Compute the per-stock moment conditions for each lag and moment type.
    Implements the inner sum over investors in eqs. (C.6) and (C.7).

    Inputs:
        base0   : ndarray (I, N), w_i * u_{i,n,t} (for lag 0, eq. C.6)
        base1   : ndarray (I, N), w_i * u_{i,n,t} * sign(Z_tilde) * z_{i,n,t} (for lags 1..L, eq. C.7)
        CB_all  : ndarray (L+1, I, N), GIV z_{i,n,t-l} for each lag
        Cts     : ndarray (L+1, I, N, P), X_{i,t,k} indicator matrices
        mf      : ndarray (L+1, I, N), data mask (bool)
        Nmom0   : int, number of zeta_1 moment conditions (lag 0)
        Nmom1   : int, number of zeta_2 moment conditions (each lag > 0)
        err_out : ndarray (L+1, max(Nmom0,Nmom1), N), output (written in place)
    """
    L, N, T = mf.shape
    for lag in nb.prange(L):
        base = base0 if lag == 0 else base1          # (N,T)
        wd   = base * CB_all[lag]                    # (N,T)

        # investor weights
        row_cnt = mf[lag].sum(axis=1)


        # Note that by not dividing by the number of holdings for this investor, we are 
        # overweighting investors with more holdings, as discussed in Appendix C.2.2.
        denom   = np.where(row_cnt > 0.0, 1.0, 1.0)   # denominator never 0

        
        
        scale   = (row_cnt > 0.0) * (1.0 / denom)         # mask zeros back in


        if lag == 0:
            start = 0
            Nmom  = Nmom0
        else:
            start = Nmom0
            Nmom  = Nmom1

        for m in range(Nmom):
            C = Cts[lag, :, :, start + m]           # (N,T)
            out = err_out[lag, m]

            # reduce over investors, parallel over stocks
            for t in range(T):
                s = 0.0
                for i in range(N):
                    s += scale[i] * wd[i, t] * C[i, t]
                out[t] = s

def get_error_mat_fast(zeta, qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list):
    """
    Compute per-stock GMM moment conditions for all lags and moment types.

    Implements eqs. (C.6) and (C.7)/(C.12) from Appendix C.2.

    Inputs:
        zeta         : ndarray (P,), current parameter vector [zeta_1; zeta_2]
        qmat         : ndarray (L+1, I, N), residualized Delta_q
        Cpts         : ndarray (L+1, I, N, P), residualized Delta_p stacked by moment type
        Cts          : ndarray (L+1, I, N, P), X_{i,t} indicator matrices
        Smat         : ndarray (L+1, I, N), lagged ownership shares
        exclude_mask : ndarray (L+1, I, N), True where data is missing
        Nmom_list    : [int, int], number of moment conditions for [lag 0, each lag > 0]

    Output:
        err : ndarray (L+1, max(Nmom_list), N), per-stock moment conditions
    """
    # Step 1: Dimensions and masks
    L, N, T = qmat.shape                                                # L+1 lags, I investors, N stocks
    mask = ~exclude_mask                                                # (L+1, I, N)
    mf = mask.astype(qmat.dtype)                                        # (L+1, I, N)

    # Step 2: Construct u_{i,n,t}(zeta) = check_q + zeta' * check_p  (eq. C.3)
    umat, Sm = build_umat_sm64(Cpts, qmat, Smat, zeta, mf)            # each (L+1, I, N)

    # Step 3: Compute precision weights w_i (eq. C.8)
    # w_i propto 1 / Var^CX(u_{i,n,t}), normalized to sum to 1
    sumsq = (umat * umat * mf).sum(axis=2)                             # (L+1, I)
    cnt = mask.sum(axis=2)                                              # (L+1, I)

    sigma_u2 = np.zeros_like(sumsq)                                     # (L+1, I)
    np.divide(sumsq, cnt, out=sigma_u2, where=(cnt > 0))
    sigma_u2 = winsorize_np(sigma_u2, limits=(0.01, 0.99))
    
    precision = np.zeros_like(sigma_u2)                                 # (L+1, I)
    np.divide(1.0, sigma_u2, out=precision, where=(sigma_u2 > 0))

    sp = precision.sum()
    if sp > 0:
        precision /= sp
    else:
        precision.fill(1.0 / N)

    # Step 4: Compute GIV z_{i,n,t-l} for each lag (leave-one-type-out, eq. 15)
    CB_minus_diagB_cache = build_CB_cache_numba(umat, Sm, block_id, use_diff=True)

    current_period_giv = CB_minus_diagB_cache[0]                        # (I, N) = z_{i,n,t}
     
    # Step 5: Compute sign(sum_l z_{i,n,t-l}) for the interaction moment (eq. C.12)
    sign_tensor = np.sign(sum(CB_minus_diagB_cache[1:]) if L > 1 else np.zeros((N, T)))  # (I, N)
    
    # Step 6: Build weighted base arrays for the two sets of moment conditions
    # base0: w_i * u_{i,n,t}                        (for eq. C.6, lag 0)
    # base1: w_i * u_{i,n,t} * sign * z_{i,n,t}     (for eq. C.7/C.12, lags > 0)
    base_arrays = [ umat[0] *  precision[0, :, None],                   # (I, N)
                    umat[0] * sign_tensor * precision[0, :, None] * current_period_giv ]  # (I, N)

    # Step 7: Compute per-stock moment conditions
    CB_all = np.stack(CB_minus_diagB_cache)                             # (L+1, I, N)

    # Demean lagged GIVs: Z_tilde = z_{i,n,t-l} - E^CX[z_{i,n,t-l}]  (eq. C.12)
    if CB_all.shape[0] <= 1:
        mean_abs_lag_GIV = np.zeros(CB_all.shape[1], dtype=CB_all.dtype)
    else:
        lag_sum = CB_all[1:].sum(axis=0)                                # (I, N)
        S = np.abs(lag_sum) / (CB_all.shape[0] - 1)

        mask0 = mf[0].astype(bool)                                     # (I, N)
        counts = mask0.sum(axis=1)                                      # (I,)
        safe_counts = np.where(counts > 0, counts, 1)
        mean_abs_lag_GIV = (S * mask0).sum(axis=1) / safe_counts       # (I,)
        mean_abs_lag_GIV = np.where(counts > 0, mean_abs_lag_GIV, 0.0)

    CB_all[1:] -= (sign_tensor * mean_abs_lag_GIV[:, None])[None, :, :]

    Nmom0, Nmom1 = Nmom_list
    err = np.empty((L, max(Nmom_list), T), dtype=base_arrays[0].dtype)  # (L+1, max(P1,P2), N)

    for arr in (base_arrays[0], base_arrays[1]):
        np.nan_to_num(arr, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        
    step7_numba(base_arrays[0], base_arrays[1],
                CB_all, Cts, mf,
                np.int64(Nmom0), np.int64(Nmom1),
                err)

    return err


def moment_conditions_fast(zeta, qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list):
    """
    Exactly-identified GMM objective: eqs. (C.6) and (C.7) averaged across stocks.

    Aggregates lagged moment conditions by summing across lags 1..L.
    Returns a vector of length (Nmom_list[0] + Nmom_list[1]) that should be
    zero at the true parameter values.

    Inputs:  same as get_error_mat_fast
    Output:  ndarray (Nmom_list[0] + Nmom_list[1],), mean moment conditions
    """
    L, N, T = qmat.shape

    err = get_error_mat_fast(zeta, qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list)

    # Stack: lag-0 moments (eq. C.6) + sum of lagged moments (eq. C.7)
    err_scaled = [
        err[0, :Nmom_list[0], :],                      # (Nmom0, N) current-period
        err[1:, :Nmom_list[1], :].sum(0)                # (Nmom1, N) sum of lagged
    ]

    err_mat = np.vstack(err_scaled)                     # (Nmom0+Nmom1, N)
    w = np.ones(err_mat.shape[1]) / err_mat.shape[1]
 
    return np.average(err_mat, axis=1, weights=w)       # (Nmom0+Nmom1,)


def moment_conditions_fast_overIdentified(zeta, qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list):
    """
    Over-identified GMM objective: eqs. (C.11) and (C.12) averaged across stocks.

    Instead of summing lagged moments across lags, keeps each lag separate.
    This gives (Nmom0 + L * Nmom1) moment conditions for P = Nmom0 + Nmom1 parameters,
    yielding overidentification that is used for root selection via the J-statistic
    (Appendix C.2.2).

    Inputs:  same as get_error_mat_fast
    Output:  ndarray (Nmom0 + L*Nmom1,), mean moment conditions
    """
    L, N, T = qmat.shape

    err = get_error_mat_fast(zeta, qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list)

    err_scaled = [
        err[0, :Nmom_list[0], :] ,
    ] \
    + \
    [
        err[ i, :Nmom_list[1], :]
        for i in range( 1, L)
    ]

    err_mat = np.vstack(err_scaled)                     # (Nmom0 + (L-1)*Nmom1, N)
    w = np.ones(err_mat.shape[1]) / err_mat.shape[1]

    return np.average(err_mat, axis=1, weights=w)       # (Nmom0 + (L-1)*Nmom1,)

def get_J_stat(zeta, qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list, obj_fun):
    """
    Compute J-statistic = sum of squared moment conditions (eq. C.13).
    Used for root selection among multiple GMM solutions (Appendix C.2.2).
    """
    f = partial(obj_fun,
                qmat=qmat, Cpts=Cpts, Cts=Cts,
                Smat=Smat, exclude_mask=exclude_mask,
                Nmom_list=Nmom_list)

    return sum(f(zeta) ** 2)

 

####################################################################################################
# STARTING POINT SELECTION FOR NONLINEAR GMM
####################################################################################################

def best_starting_points(obj_fun, bounds, args, N=10, oversample=50_000):
    """
    Find good starting points for nonlinear GMM using Sobol quasi-random sampling.

    The GMM objective (C.10) can have multiple local minima/roots. This function
    evaluates the objective at many quasi-random candidate points and returns
    the N candidates with the smallest sum of squared moment conditions.

    Inputs:
        obj_fun    : callable, moment condition function (returns residual vector)
        bounds     : (lo, hi), each ndarray (P,), parameter search bounds
        args       : tuple, additional arguments passed to obj_fun
        N          : int, number of best starting points to return
        oversample : int, total number of candidates to evaluate

    Outputs:
        X[best_idx] : ndarray (N, P), best starting points
        X           : ndarray (oversample, P), all candidates
        residuals   : ndarray, last evaluated residual vector
        losses      : ndarray (oversample,), sum of squared residuals at each candidate
    """
    # Step 1: Generate Sobol quasi-random candidates in the parameter hypercube
    lo, hi = map(np.asarray, bounds)
    d = lo.size                                                         # P = number of parameters
    
    rng = qmc.Sobol(d, scramble=True, seed=0)
    np.random.default_rng(0)
    
    X_unit = rng.random_base2(int(np.ceil(np.log2(oversample))))        # (>=oversample, P) in [0,1]^P
    X = lo + (hi - lo) * X_unit[:oversample]                            # (oversample, P)
    
    # Step 2: Evaluate ||moment_conditions(zeta)||^2 at each candidate
    losses = np.empty(oversample)                                       # (oversample,)
    t0 = time.time()
    
    for i, x in enumerate(X):
        residuals = obj_fun(x, *args)
        losses[i] = np.dot(residuals, residuals)
    
    # Step 3: Return N candidates with smallest loss
    best_idx = np.argpartition(losses, N)[:N]
    
    return X[best_idx], X,  residuals, losses


####################################################################################
# GMM ESTIMATION: SOLVE FOR ZETA = [zeta_{1,0}, zeta_{1,AS}, zeta_{2,0}, zeta_{2,AS}]
####################################################################################
# Minimize the exactly-identified GMM objective (C.10):
#   min_{zeta} sum_k ( sum_i w_i * E^CX[X_{i,k} * u * z] )^2
#            + sum_k ( sum_i w_i * E^CX[X_{i,k} * u * Z_tilde * z] )^2
#
# We use multiple starting points (Sobol-sampled) and solve with fsolve.
# The best root is selected using overidentifying restrictions (Appendix C.2.2).
####################################################################################

print('Starting GMM estimation optimization...')
t0 = time.time()

# Parameter vector: [zeta_{1,0}, zeta_{1,ActiveShare}, zeta_{2,0}, zeta_{2,ActiveShare}]
initial_guess = np.array([.5] * Nmom_list[0] + [0.00] * Nmom_list[1])  # (P,)

# Search bounds for starting point generation (not imposed on fsolve itself)
bounds = (
    np.array([.01] * Nmom_list[0] + [-.1] * Nmom_list[1]),
    np.array([1] * Nmom_list[0] + [.1] * Nmom_list[1])
)

N = 20

t0 = time.time()
starts, all_vecs, residuals, losses  = best_starting_points(
    moment_conditions_fast,
    bounds=bounds,
    args=(qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list),
    N=N,
    oversample=500
)

sorted_indices = np.argsort(losses)[:N]
starts = all_vecs[sorted_indices]                                       # (N, P)
starts = np.vstack([starts, np.clip(initial_guess, bounds[0], bounds[1])])  # (N+1, P)

min_cost = 10
min_sse = 10
best_start_guess_index = 0
l = []
l_fsolve = []
l_ls = []




def process_single_start(args):
    """
    Solve the GMM system from one starting point using fsolve (root finding).

    Also computes J-statistics under both the exactly-identified and
    over-identified moment conditions for root selection (Appendix C.2.2).

    Inputs:
        args : tuple containing:
            start_guess_index : int
            start_guess       : ndarray (P,)
            qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list : data tensors

    Output:
        dict with keys:
            'start_guess_index' : int
            'fsolve_solution'   : tuple from scipy.optimize.fsolve (full_output=True)
            'fsolve_sse'        : float, sum of squared moment conditions at solution
    """
    (start_guess_index, start_guess, qmat, Cpts, Cts, Smat, exclude_mask, 
    Nmom_list) = args
    
    t0 = time.time() 
    this_solution = fsolve(
        moment_conditions_fast,
        start_guess,
        args=(qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list),
        full_output=True,
        epsfcn=1e-10, 
    )

    # Compute J-statistics for root selection
    this_solution[1]['J_stats'] = {}
    this_solution[1]['J_stats']['moment_conditions_fast'] = get_J_stat(this_solution[0], qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list, moment_conditions_fast)
    this_solution[1]['J_stats']['moment_conditions_fast_overIdentified'] = get_J_stat(this_solution[0], qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list, moment_conditions_fast_overIdentified)

    mask =  deepcopy(qmat)
    mask[:] = 1
    for data_mat in [qmat, Cpts[:, :, :, 0], Cts[:, :, :, 0], Smat]:
        mask[np.isnan(data_mat)] *= 0
        mask[ data_mat == 0 ] *= 0

    this_solution[1]['nobs']  =  mask.sum()
    this_solution[1]['qmat_shape'] = qmat.shape

    fsolve_sse = sum(this_solution[1]['fvec'] ** 2)

    return {
        'start_guess_index': start_guess_index,
        'fsolve_solution': this_solution,
        'fsolve_sse': fsolve_sse,
    }
 
####################################################################################
# PARALLEL EXECUTION: SOLVE FROM MULTIPLE STARTING POINTS
####################################################################################

args_list = [
    (i, starts[i], qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list)
    for i in range(len(starts))
]

l_fsolve = []
l_ls = []
min_sse = float('inf')
min_cost = float('inf')
solution = None
result = None
best_start_guess_index = None

num_workers = 8

num_workers = 1

if not local:
    num_workers = int(os.environ['SLURM_JOB_CPUS_PER_NODE']) - 1
    num_workers = min(num_workers, 6)
    num_workers = 1

output_folder = output_path[ : output_path.rfind('/') + 1]
os.makedirs(output_folder, exist_ok=True)

def save_results():
    """
    Checkpoint current best results to disk as a pickle file.
    Called after each starting point completes for fault tolerance.
    """

    try:
        solution is None
    
    except:
        best_converged_param = -float('inf')
        converged_solution = None
        
        for fsolve_result in l_fsolve:
            if fsolve_result is not None and len(fsolve_result) > 1:
                if  fsolve_result[2] == 1:
                    if fsolve_result[0][0] > best_converged_param:
                        best_converged_param = fsolve_result[0][0]
                        converged_solution = fsolve_result
        
        if converged_solution is not None:
            solution = converged_solution
            
        else:
            best_sse = float('inf')
            for fsolve_result in l_fsolve:

                    fsolve_sse = sum(fsolve_result[1]['fvec'] ** 2)
                    if fsolve_sse < best_sse:
                        best_sse = fsolve_sse
                        solution = fsolve_result
    
    
    with open(output_folder + str(output_quarters[-1]) + '.p', 'wb') as f:
        pickle.dump([solution, result, l_fsolve, l_ls, starts], f)

# Run fsolve from each starting point, collect results as they complete
with ThreadPoolExecutor(max_workers=num_workers) as executor:
    future_to_index = {
        executor.submit(process_single_start, args): args[0] 
        for args in args_list
    }
    
    found_convergence = False
    
    for future in concurrent.futures.as_completed(future_to_index):
        res = future.result()
        l_fsolve.append(res['fsolve_solution'])
        
        updated = False
        
        # Keep best converged root with zeta_{1,0} > 0 (demand curves slope down)
        if res['fsolve_sse'] < min_sse and res['fsolve_solution'][0][0] > 0 and res['fsolve_solution'][2] == 1:
            min_sse = res['fsolve_sse']
            solution = res['fsolve_solution']
            updated = True
        
        save_results()
        completed = len(l_fsolve)
        total = len(args_list)

save_results() 

print("=" * 60)
print("FINAL OPTIMIZATION RESULTS")
print("=" * 60)
print("Best fsolve solution:")
print(solution)
    
 
####################################################################################
# GIV CALCULATION AND OUTPUT
####################################################################################
# Using the estimated zeta, construct the GIV z_{i,n,t} = sum_{j!=i} S_{j,n,t-1} * u_{j,n,t}
# and output investor-stock-quarter level GIV values as a .dta file.
####################################################################################

print('Starting GIV generation...')

output_folder = output_path[ : output_path.rfind('/') + 1]

exclude_mask = np.isnan(qmat) + np.isnan(Smat)
Smat[np.isnan(Smat)] = 0.0
mask = ~exclude_mask
mf = mask.astype(qmat.dtype)

# Reconstruct leave-one-type-out matrix
type_pd = pd.DataFrame(type_tensor[0])
type_vec = type_pd.T.ffill().T.iloc[:, -1]
vals = type_vec.to_numpy()

type_neq_mat = (vals[:, None] != vals).astype(int)                      # (I, I)
neq_mat = type_neq_mat 


def get_giv(zeta, qmat, Cpts, Smat, exclude_mask):
    """
    Compute the GIV z_{i,n,t} for each lag using estimated zeta.

    z_{i,n,t-l} = sum_{j: type(j) != type(i)} S_{j,n,t-1-l} * u_{j,n,t-l}(zeta)

    where u_{j,n,t} = check_q_{j,n,t} + zeta' * check_p_{n,t}

    Inputs:
        zeta         : ndarray (P,), estimated elasticity parameters
        qmat         : ndarray (L+1, I, N), residualized Delta_q
        Cpts         : ndarray (L+1, I, N, P), residualized Delta_p by moment type
        Smat         : ndarray (L+1, I, N), lagged ownership shares
        exclude_mask : ndarray (L+1, I, N), True where data is missing

    Outputs:
        giv_matrices : list of length (L+1), each ndarray (I, N) = z_{i,n,t-l}
        precision    : ndarray (L+1, I), precision weights w_i for each lag
    """
    # Step 1: Dimensions and masks
    L, N, T = qmat.shape
    num_zeta_elems = zeta.shape[0]
    Nmom = zeta.shape[0] // 2
    
    mask = ~exclude_mask                                                # (L+1, I, N)
    mf = mask.astype(qmat.dtype)                                        # (L+1, I, N)

    # Step 2: Construct u_{i,n,t}(zeta) = check_q + zeta' * check_p
    umat = qmat * mf                                                   # (L+1, I, N)
    for m in range(num_zeta_elems): 
        umat += Cpts[:, :, :, m] * zeta[m] * mf
    umat = np.nan_to_num(umat, nan=0.0) * mf

    # Step 3: Clean ownership shares
    Sm = np.nan_to_num(Smat, nan=0.0) * mf                             # (L+1, I, N)

    # Step 4: Compute precision weights (eq. C.8)
    sumsq = (umat * umat * mf).sum(axis=2)                             # (L+1, I)
    cnt = mask.sum(axis=2).astype(qmat.dtype)                           # (L+1, I)
    sigma_u2 = np.where(cnt > 0, sumsq / cnt, 0.0)
    sigma_u2 = winsorize_np(sigma_u2, limits=(0.01, 0.99))
    
    precision = np.where(sigma_u2 > 0, 1.0 / sigma_u2, 0.0)           # (L+1, I)
    sp = precision.sum()
    if sp > 0:
        precision /= sp
    else:
        precision.fill(1.0 / N)

    # Step 5: Compute GIV for each lag: z_{i,n,t-l} = neq_mat @ (u * S)
    giv_matrices = []
    for lag in range(0, L): 
        B = umat[lag] * Sm[lag]                                         # (I, N)
        CB = neq_mat @ B                                                # (I, N)
        giv_matrices.append(CB)

    return giv_matrices, precision



def get_CB_all(zeta, qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list):
    """
    Compute per-stock GMM moment conditions for all lags and moment types.

    Implements eqs. (C.6) and (C.7)/(C.12) from Appendix C.2.

    Inputs:
        zeta         : ndarray (P,), current parameter vector [zeta_1; zeta_2]
        qmat         : ndarray (L+1, I, N), residualized Delta_q
        Cpts         : ndarray (L+1, I, N, P), residualized Delta_p stacked by moment type
        Cts          : ndarray (L+1, I, N, P), X_{i,t} indicator matrices
        Smat         : ndarray (L+1, I, N), lagged ownership shares
        exclude_mask : ndarray (L+1, I, N), True where data is missing
        Nmom_list    : [int, int], number of moment conditions for [lag 0, each lag > 0]

    Output:
        err : ndarray (L+1, max(Nmom_list), N), per-stock moment conditions
    """
    # Step 1: Dimensions and masks
    L, N, T = qmat.shape                                                # L+1 lags, I investors, N stocks
    mask = ~exclude_mask                                                # (L+1, I, N)
    mf = mask.astype(qmat.dtype)                                        # (L+1, I, N)

    # Step 2: Construct u_{i,n,t}(zeta) = check_q + zeta' * check_p  (eq. C.3)
    umat, Sm = build_umat_sm64(Cpts, qmat, Smat, zeta, mf)            # each (L+1, I, N)

    # Step 3: Compute precision weights w_i (eq. C.8)
    # w_i propto 1 / Var^CX(u_{i,n,t}), normalized to sum to 1
    sumsq = (umat * umat * mf).sum(axis=2)                             # (L+1, I)
    cnt = mask.sum(axis=2)                                              # (L+1, I)

    sigma_u2 = np.zeros_like(sumsq)                                     # (L+1, I)
    np.divide(sumsq, cnt, out=sigma_u2, where=(cnt > 0))
    sigma_u2 = winsorize_np(sigma_u2, limits=(0.01, 0.99))
    
    precision = np.zeros_like(sigma_u2)                                 # (L+1, I)
    np.divide(1.0, sigma_u2, out=precision, where=(sigma_u2 > 0))

    sp = precision.sum()
    if sp > 0:
        precision /= sp
    else:
        precision.fill(1.0 / N)

    # Step 4: Compute GIV z_{i,n,t-l} for each lag (leave-one-type-out, eq. 15)
    CB_minus_diagB_cache = build_CB_cache_numba(umat, Sm, block_id, use_diff=True)

    current_period_giv = CB_minus_diagB_cache[0]                        # (I, N) = z_{i,n,t}
     
    # Step 5: Compute sign(sum_l z_{i,n,t-l}) for the interaction moment (eq. C.12)
    sign_tensor = np.sign(sum(CB_minus_diagB_cache[1:]) if L > 1 else np.zeros((N, T)))  # (I, N)
    
    # Step 6: Build weighted base arrays for the two sets of moment conditions
    # base0: w_i * u_{i,n,t}                        (for eq. C.6, lag 0)
    # base1: w_i * u_{i,n,t} * sign * z_{i,n,t}     (for eq. C.7/C.12, lags > 0)
    base_arrays = [ umat[0] *  precision[0, :, None],                   # (I, N)
                    umat[0] * sign_tensor * precision[0, :, None] * current_period_giv ]  # (I, N)

    # Step 7: Compute per-stock moment conditions
    CB_all = np.stack(CB_minus_diagB_cache)                             # (L+1, I, N)

    # Demean lagged GIVs: Z_tilde = z_{i,n,t-l} - E^CX[z_{i,n,t-l}]  (eq. C.12)
    if CB_all.shape[0] <= 1:
        mean_abs_lag_GIV = np.zeros(CB_all.shape[1], dtype=CB_all.dtype)
    else:
        lag_sum = CB_all[1:].sum(axis=0)                                # (I, N)
        S = np.abs(lag_sum) / (CB_all.shape[0] - 1)

        mask0 = mf[0].astype(bool)                                     # (I, N)
        counts = mask0.sum(axis=1)                                      # (I,)
        safe_counts = np.where(counts > 0, counts, 1)
        mean_abs_lag_GIV = (S * mask0).sum(axis=1) / safe_counts       # (I,)
        mean_abs_lag_GIV = np.where(counts > 0, mean_abs_lag_GIV, 0.0)

    CB_all[1:] -= (sign_tensor * mean_abs_lag_GIV[:, None])[None, :, :]
    
    
    return CB_all, precision, sign_tensor



# Load estimation results and select the best root using J-statistic (Appendix C.2.2):
# 1. Among converged roots with zeta_{1,0} > 0, pick the one minimizing
#    the over-identified J-statistic (eq. C.13).
# 2. If no converged root has zeta_{1,0} > 0, fall back to all converged roots.
# 3. If no root converged, then use the root with the lowest J-statistic.
all_results = pickle.load(open(output_folder + str(output_quarters[-1]) + '.p', 'rb'))


all_converged_roots =  sorted([ (s[1]['J_stats']['moment_conditions_fast_overIdentified'], s)  for s in all_results[2] if s[2] == 1 and s[0][0] > 0], key = lambda s : s[0])

if len(all_converged_roots) <= 0:
    all_converged_roots =  sorted([ (s[1]['J_stats']['moment_conditions_fast_overIdentified'], s)  for s in all_results[2] if s[2] == 1], key = lambda s : s[0])

if len(all_converged_roots) <= 0:
    all_converged_roots =  sorted([ (s[1]['J_stats']['moment_conditions_fast'], s)  for s in all_results[2] ], key = lambda s : s[0])

result_fsolve = all_converged_roots[0][1]
zeta = result_fsolve[0]                                             # (P,) = [zeta_{1,0}, zeta_{1,AS}, zeta_{2,0}, zeta_{2,AS}]

# Compute GIV at the selected zeta
# list_of_matrices, precision = get_giv(zeta, qmat, Cpts, Smat, exclude_mask)

list_of_matrices, precision, sign_tensor = get_CB_all(zeta, qmat, Cpts, Cts, Smat, exclude_mask, Nmom_list)

# Convert GIV matrices from (I x N) arrays back to long-format DataFrames
# (mgrno, permno, quarter) with GIV values for each lag.
list_of_giv_dfs = []
list_of_delta_q_dfs = []
list_of_type_dfs = []
list_of_Cpts_dfs = []
list_of_Cpts_interaction_dfs = []
list_of_Cts_dfs = []

for lag in range(num_lagged_returns + 1):
    this_lag_giv = list_of_matrices[lag]                                # (I, N) GIV for this lag
    this_lag_giv[mf[lag] == 0] = np.nan

    # Reconstruct column (stock) labels from the concatenated quarterly DataFrames
    column_reference_list = []
    
    for qtr in range(num_quaters):
        column_reference_list.append(qmat_pd_list[qtr][lag])

    mat = pd.concat(column_reference_list, axis=1)
    this_lag_giv_df = pd.DataFrame(this_lag_giv, index=mat.index)       # (I, N) DataFrame

    # Split back into per-quarter long-format DataFrames
    quarterly_giv_list = []
    quarterly_delta_q_list = []
    quarterly_type_list = []
    quarterly_Cpts_list = []
    quarterly_Cpts_interaction_list = []
    quarterly_Cts_list = []

    start_col = 0

    for qtr in range(num_quaters):
        num_stocks_this_qtr = len(qmat_pd_list[qtr][lag].columns)

        this_lag_this_quarter_delta_q_df = pd.melt(qmat_pd_list[qtr][lag], ignore_index=False).reset_index()
        this_lag_this_quarter_delta_q_df = this_lag_this_quarter_delta_q_df.rename( columns = {'value' : 'delta_q_lag_%s' % lag})
        this_lag_this_quarter_delta_q_df['qtr'] = job_file_index + qtr
        quarterly_delta_q_list.append(this_lag_this_quarter_delta_q_df)

        this_lag_this_quarter_type_df = pd.melt(type_pd_list[qtr][lag], ignore_index=False).reset_index()
        this_lag_this_quarter_type_df = this_lag_this_quarter_type_df.rename( columns = {'value' : 'type_lag_%s' % lag})
        this_lag_this_quarter_type_df['qtr'] = job_file_index + qtr
        quarterly_type_list.append(this_lag_this_quarter_type_df)

        # price_mat = np.sum(Cpts_pd_list[qtr][lag], 0)
        this_quarter_list_of_Cpts_dfs = []
        for index in range(len(Cpts_pd_list[qtr][lag])):

            price_mat = Cpts_pd_list[qtr][lag][index]
                
            price_mat = pd.DataFrame(price_mat, index=Cpts_pd_list[qtr][lag][0].index, columns=Cpts_pd_list[qtr][lag][0].columns)

            this_lag_this_quarter_Cpts_df = pd.melt(price_mat, ignore_index=False).reset_index()
            if index == 0:
                this_lag_this_quarter_Cpts_df = this_lag_this_quarter_Cpts_df.rename( columns = {'value' : 'Cpts_lag_%s' % lag})
            else:
                this_lag_this_quarter_Cpts_df = this_lag_this_quarter_Cpts_df.rename( columns = {'value' : 'Cpts_lag_%s_x_AS' % lag})
            this_lag_this_quarter_Cpts_df['qtr'] = job_file_index + qtr

            this_lag_this_quarter_Cpts_df = this_lag_this_quarter_Cpts_df.set_index(['mgrno', 'permno', 'qtr'])

            this_quarter_list_of_Cpts_dfs.append(this_lag_this_quarter_Cpts_df)

        
        this_lag_this_quarter_Cpts_df = pd.concat(this_quarter_list_of_Cpts_dfs, axis=1)
        quarterly_Cpts_list.append(this_lag_this_quarter_Cpts_df)

        # price_mat = np.sum(Cpts_pd_list[qtr][lag], 0)
        this_quarter_list_of_Cts_dfs = []
        for index in range(len(Cts_pd_list[qtr][lag])):

            this_Cts_mat = Cts_pd_list[qtr][lag][index]
                
            this_Cts_mat = pd.DataFrame(this_Cts_mat, index=Cts_pd_list[qtr][lag][0].index, columns=Cts_pd_list[qtr][lag][0].columns)

            this_lag_this_quarter_Cts_df = pd.melt(this_Cts_mat, ignore_index=False).reset_index()
            if index == 0:
                this_lag_this_quarter_Cts_df = this_lag_this_quarter_Cts_df.rename( columns = {'value' : 'Cts_lag_%s' % lag})
            else:
                this_lag_this_quarter_Cts_df = this_lag_this_quarter_Cts_df.rename( columns = {'value' : 'Cts_lag_%s_x_AS' % lag})
            this_lag_this_quarter_Cts_df['qtr'] = job_file_index + qtr

            this_lag_this_quarter_Cts_df = this_lag_this_quarter_Cts_df.set_index(['mgrno', 'permno', 'qtr'])

            this_quarter_list_of_Cts_dfs.append(this_lag_this_quarter_Cts_df)

        
        this_lag_this_quarter_Cts_df = pd.concat(this_quarter_list_of_Cts_dfs, axis=1)
        quarterly_Cts_list.append(this_lag_this_quarter_Cts_df)


        # this_lag_this_quarter_Cts_df = pd.melt(Cts_pd_list[qtr][lag], ignore_index=False).reset_index()
        # this_lag_this_quarter_Cts_df = this_lag_this_quarter_Cts_df.rename( columns = {'value' : 'Cts_lag_%s' % lag})
        # this_lag_this_quarter_Cts_df['qtr'] = job_file_index + qtr
        # quarterly_Cts_list.append(this_lag_this_quarter_Cts_df)

        
        # int_price_mat = np.sum(Cpts_interaction_pd_list[qtr][lag], 0)
        this_quarter_list_of_Cpts_interaction_dfs = []
        for index in range(len(Cpts_interaction_pd_list[qtr][lag])): 
            int_price_mat = Cpts_interaction_pd_list[qtr][lag][index]
            int_price_mat = pd.DataFrame(int_price_mat, index=Cpts_interaction_pd_list[qtr][lag][0].index, columns=Cpts_interaction_pd_list[qtr][lag][0].columns)

            this_lag_this_quarter_Cpts_interaction_df = pd.melt(int_price_mat, ignore_index=False).reset_index()
            if index == 0:
                this_lag_this_quarter_Cpts_interaction_df = this_lag_this_quarter_Cpts_interaction_df.rename( columns = {'value' : 'Cpts_interaction_lag_%s' % lag})
            else:
                this_lag_this_quarter_Cpts_interaction_df = this_lag_this_quarter_Cpts_interaction_df.rename( columns = {'value' : 'Cpts_interaction_lag_%s_x_AS' % lag})
            this_lag_this_quarter_Cpts_interaction_df['qtr'] = job_file_index + qtr

            this_lag_this_quarter_Cpts_interaction_df = this_lag_this_quarter_Cpts_interaction_df.set_index(['mgrno', 'permno', 'qtr'])

            this_quarter_list_of_Cpts_interaction_dfs.append(this_lag_this_quarter_Cpts_interaction_df)
        
        this_lag_this_quarter_Cpts_interaction_df = pd.concat(this_quarter_list_of_Cpts_interaction_dfs, axis=1)
        quarterly_Cpts_interaction_list.append(this_lag_this_quarter_Cpts_interaction_df)
            
        this_lag_this_quarter_giv_df = this_lag_giv_df.loc[
            qmat_pd_list[qtr][lag].index
        ].iloc[:, start_col:start_col + num_stocks_this_qtr]

        this_lag_this_quarter_giv_df.columns = list(mat.columns)[
            start_col:start_col + num_stocks_this_qtr
        ]

        start_col_old = start_col
        start_col = start_col + num_stocks_this_qtr

        # Melt from wide (I x N_t) to long (observations x 1)
        this_lag_this_quarter_giv_df = pd.melt(
            this_lag_this_quarter_giv_df, 
            ignore_index=False
        ).reset_index()
        this_lag_this_quarter_giv_df.columns = ['mgrno', 'permno', 'GIV_lag_%s' % lag]
        this_lag_this_quarter_giv_df['qtr'] = job_file_index + qtr


        if lag == 0:
            CB_all = np.stack(list_of_matrices)

            Z_gmm_interaction = np.zeros_like(list_of_matrices[0])
            Z_gmm_interaction_AS = np.zeros_like(list_of_matrices[0])

            for l in range(1, num_lagged_returns + 1):
                Z_gmm_interaction += (
                    list_of_matrices[0]
                    * sign_tensor
                    * CB_all[l]
                    * Cts[l, :, :, Nmom_list[0] + 0]
                )

                Z_gmm_interaction_AS += (
                    list_of_matrices[0]
                    * sign_tensor
                    * CB_all[l]
                    * Cts[l, :, :, Nmom_list[0] + 1]
                )

            Z_level = list_of_matrices[0] * Cts[0, :, :, 0]
            Z_level_AS = list_of_matrices[0] * Cts[0, :, :, 1]

            for Zmat, zname in [
                (Z_level, "GIV_gmm_level"),
                (Z_level_AS, "GIV_gmm_level_x_AS"),
                (Z_gmm_interaction, "GIV_gmm_interaction"),
                (Z_gmm_interaction_AS, "GIV_gmm_interaction_x_AS"),
            ]:
                Z_q = pd.DataFrame(
                    Zmat[:, start_col_old:start_col],
                    index=mat.index,
                    columns=list(mat.columns)[start_col_old:start_col],
                )

                Z_q = pd.melt(Z_q, ignore_index=False).reset_index()
                Z_q.columns = ["mgrno", "permno", zname]
                Z_q["qtr"] = job_file_index + qtr

                this_lag_this_quarter_giv_df = this_lag_this_quarter_giv_df.merge(
                    Z_q,
                    on=["mgrno", "permno", "qtr"],
                    how="left",
                )

            # sys.exit(0)
        quarterly_giv_list.append(this_lag_this_quarter_giv_df)

    # Combine all quarters for this lag and set multi-index
    this_lag_combined = pd.concat(quarterly_giv_list, axis=0)
    this_lag_combined = this_lag_combined.set_index(['mgrno', 'permno', 'qtr'])

    this_lag_combined_delta_q = pd.concat(quarterly_delta_q_list, axis=0)
    this_lag_combined_delta_q = this_lag_combined_delta_q.set_index(['mgrno', 'permno', 'qtr'])

    this_lag_combined_type = pd.concat(quarterly_type_list, axis=0)
    this_lag_combined_type = this_lag_combined_type.set_index(['mgrno', 'permno', 'qtr'])

    this_lag_combined_Cpts = pd.concat(quarterly_Cpts_list, axis=0)
    # this_lag_combined_Cpts = this_lag_combined_Cpts.set_index(['mgrno', 'permno', 'qtr'])

    this_lag_combined_Cpts_interaction = pd.concat(quarterly_Cpts_interaction_list, axis=0)
    # this_lag_combined_Cpts_interaction = this_lag_combined_Cpts_interaction.set_index(['mgrno', 'permno', 'qtr'])

    this_lag_combined_Cts = pd.concat(quarterly_Cts_list, axis=0)
    # this_lag_combined_Cts = this_lag_combined_Cts.set_index(['mgrno', 'permno', 'qtr'])

    list_of_giv_dfs.append(this_lag_combined)
    list_of_delta_q_dfs.append(this_lag_combined_delta_q)
    list_of_type_dfs.append(this_lag_combined_type)
    list_of_Cpts_dfs.append(this_lag_combined_Cpts)
    list_of_Cpts_interaction_dfs.append(this_lag_combined_Cpts_interaction)
    list_of_Cts_dfs.append(this_lag_combined_Cts)
    
# Merge all lags into a single DataFrame keyed by (mgrno, permno, qtr)
full_giv_df = pd.concat(list_of_giv_dfs, axis=1)
full_giv_df = full_giv_df.sort_index()
full_giv_df = full_giv_df.rename(columns={'mgrno': 'bin'})


full_delta_q_df = pd.concat(list_of_delta_q_dfs, axis=1)
full_delta_q_df = full_delta_q_df.sort_index()
full_delta_q_df = full_delta_q_df.rename(columns={'mgrno': 'bin'})

full_type_df = pd.concat(list_of_type_dfs, axis=1)
full_type_df = full_type_df.sort_index()
full_type_df = full_type_df.rename(columns={'mgrno': 'bin'})

full_Cpts_df = pd.concat(list_of_Cpts_dfs, axis=1)
full_Cpts_df = full_Cpts_df.sort_index()
full_Cpts_df = full_Cpts_df.rename(columns={'mgrno': 'bin'})

full_Cpts_interaction_df = pd.concat(list_of_Cpts_interaction_dfs, axis=1)
full_Cpts_interaction_df = full_Cpts_interaction_df.sort_index()

full_Cts_df = pd.concat(list_of_Cts_dfs, axis=1)
full_Cts_df = full_Cts_df.sort_index()
full_Cts_df = full_Cts_df.rename(columns={'mgrno': 'bin'})

full_df = pd.concat([full_giv_df, full_delta_q_df, full_type_df, full_Cpts_df, full_Cpts_interaction_df, full_Cts_df], axis=1)
# full_df = full_df.dropna()


giv_lag_columns = ['GIV_lag_%s' % lag for lag in range(1, num_lagged_returns + 1)]
full_df['GIV_lag_sum'] = full_df[giv_lag_columns].sum(1) 

full_df['GIV_lag_sum_abs'] = np.abs(full_df['GIV_lag_sum'])
full_df['GIV_x_GIV_lag_sum_abs'] = full_df['GIV_lag_0'] * full_df['GIV_lag_sum_abs']

# Attach precision weights (inverse variance, eq. C.8) to each investor
precision_weight_df = pd.DataFrame(precision[0, :, None], index =  mat.index ).reset_index()
precision_weight_df = precision_weight_df.rename(columns = {0 : 'precision_weight', 'index' : 'mgrno'})

full_df = pd.merge(full_df.reset_index(), precision_weight_df, on = 'mgrno', how = 'left')
full_df = full_df.set_index(['mgrno', 'permno', 'qtr'])
 
needed_cols = [
    "delta_q_lag_0",
    "Cpts_lag_0",
    "Cpts_lag_0_x_AS",
    "Cpts_interaction_lag_0",
    "Cpts_interaction_lag_0_x_AS",
    "GIV_gmm_level",
    "GIV_gmm_level_x_AS",
    "GIV_gmm_interaction",
    "GIV_gmm_interaction_x_AS",
    "precision_weight",
]

full_df = full_df.dropna(subset=needed_cols)

####################################################################################
# WRITE GIV OUTPUT TO DISK
####################################################################################
scratch_dir = '/fs/scratch/PAS2771/'
if local:
    scratch_dir = ''

output_folder = giv_output_path[ : giv_output_path.rfind('/') + 1]

full_df.to_stata(output_folder + 'GIV_quarter_%s.dta' % output_quarters[-1])
print("GIV file written to: ", output_folder + 'GIV_quarter_%s.dta' % output_quarters[-1])

