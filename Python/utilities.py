"""
@author: Scott Campit
"""

import pandas as pd
import numpy as np

def convert_to_string(df, col):
    df[col] = df[col].fillna(-1)
    df[col] = df[col].astype(int)
    df[col] = df[col].astype(str)
    df[col] = df[col].replace('-1', np.nan)
    return df

def sep_object(df, col, regex):
    s = df[col].str.split(regex).apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = col
    del df[col]
    s = pd.DataFrame(s)
    df = pd.merge(df, s, how='inner', left_index=True, right_index=True)
    df = df.drop_duplicates(keep='first')
    return df
