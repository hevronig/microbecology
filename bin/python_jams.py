# Unordered list of useful data wrangling Python snippets. Mainly pandas, numpy, matplotlib, etc.

# Create a Python list with paths of files from a directory using IPython magic
files = !ls path/*.csv

# Import multiple files with pandas and map
dfs = list(map(lambda df: pd.read_csv(df, names=[df]), files))

# Merge multiple pandas data frames into one data frame (Python3: from functools import reduce
dfs_merged = reduce(lambda df1, df2: pd.merge(df1, df2, how='outer', left_index=True, right_index=True), dfs)

# Move legend box outside plot (by [ShitalShah](https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot))
df.myCol.plot().legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Plot data frame categorical data above/below some threshold
df['column_of_interes'].value_counts().loc[df['column_of_interes'].value_counts().get_values() > 5].plot(kind='bar');
