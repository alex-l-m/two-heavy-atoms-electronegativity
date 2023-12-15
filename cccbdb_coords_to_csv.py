import json
import pandas as pd
from os import environ
from os.path import join

env_var_name = 'CCCBDB_PATH'
# Retrieve database path from environment variable
database_path = environ[env_var_name]
coord_db_filename = 'srd101_ExpCartesians.json'
inpath = join(database_path, coord_db_filename)

# Read a list of data records, each representing an individual atom
database_json = json.load(open(inpath))['ExpCartesians']

# Create a dataframe with each record as a row
database_df = pd.DataFrame(database_json)

# Write as a csv file
outpath = 'srd101_ExpCartesians.csv.gz'
database_df.to_csv(outpath, index=False, compression='gzip')
