# Here I am using Python as it has better data handling capabilities for the data format of ERA5 output
import xarray as xr
import pandas as pd
import numpy as np
import cfgrib

# Import packages for conversting pandas dataframe to RDS format
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# Open ERA5 dataset as GRIB
ds = cfgrib.open_datasets('ERA5/adaptor.mars.internal-1717355356.2162263-17481-13-2d67889f-d258-46d8-8806-001b74a90fbb.grib')

# Transform and extract data and specific coordinate
ds_0 = ds[0].to_array().to_pandas().T
# select value at lon -64.8 due to land-sea mask parameter
ds_1 = ds[1].sel(longitude=-64.8).to_array().to_pandas().T
ds_2 = ds[2].sel(longitude=-64.8).to_array().to_pandas().T

# Match the indices
ds_1.index = ds_0.index
ds_2.index = ds_0.index

# Combine data
ds_merged = pd.concat([ds_0,ds_1,ds_2], axis=1)
ds_merged_2 = ds_merged.reset_index()

# source: https://stackoverflow.com/a/62947565


# Convert pandas dataframe to R dataframe
with localconverter(robjects.default_converter + pandas2ri.converter):
    r_df = robjects.conversion.py2rpy(ds_merged_2)

# Save R dataframe as .rds file
r_file = "processed/ERA5_data.RDS"
robjects.r.assign("ERA5", r_df)
robjects.r(f"saveRDS(ERA5, '{r_file}')")
