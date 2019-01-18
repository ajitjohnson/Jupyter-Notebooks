# Preprocess the data: Normalization
# Author: Ajit Johnson Nirmal
# Last updated: 1/15/19

import numpy as np
import pandas as pd
def mi_pp_normalize (data, qq = 0.001):
    print ("Applying log transformation and quantile normalization...")
    # log10-transform of users data
    #log_data = np.log10(data)
    log_data = data.apply(np.log10, axis=0)
    log_data = log_data.replace([np.inf, -np.inf], np.nan).fillna(0)
    # Quantile normalize of each channel of users data
    lo = np.quantile(log_data, qq, axis = 0)
    hi = np.quantile(log_data, 1 - qq, axis = 0)
    normalized_data = 2 * (log_data - lo) / (hi - lo) - 1
    return(normalized_data)
