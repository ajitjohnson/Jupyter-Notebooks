# Intial formating of data from segmentation to data analysis
# Author: Ajit Johnson Nirmal
# Last updated: 1/17/19

def mi_pp_cleandata (data):
    # Identify batches
    data_median = data.filter(regex='Median')
    # Drop unwanted DAPI channels and background chennels
    data_median = data_median.drop(list(data_median.filter(regex='DNA |back', axis=1)), axis=1)
    # Change column name
    data_median.columns = data_median.columns.str.split('_').str[0]
    return (data_median)
