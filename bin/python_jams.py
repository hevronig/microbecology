# Unordered list of useful data wrangling Python snippets. Mainly pandas, numpy, matplotlib, etc.

# Create a Python list with paths of files from a directory using IPython magic
files = !ls path/*.csv 

# Import multiple files with pandas and map
dfs = list(map(lambda df: pd.read_csv(df, names=[df]), files))

# Merge multiple pandas data frames into one data frame (Python3: from functools import reduce
dfs_merged = reduce(lambda df1, df2: pd.merge(df1, df2, how='outer', left_index=True, right_index=True), dfs)

#


#
