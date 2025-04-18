# -*- coding: utf-8 -*-

import pandas as pd
import glob
import os

directory = "/path/to/directory"

# Get a list of all .txt files in the directory
txt_files = glob.glob("*.txt")

# Create an empty dictionary to store the DataFrames
dfs = {}

# Loop through each .txt file
for file in txt_files:
    # Extract the DataFrame name from the file name
    df_name = os.path.basename(file).split('_reads_counts.txt')[0]

    # Read the file into a DataFrame
    df = pd.read_csv(file, delimiter='\t', header=None, names=['Data'])

    # Split the first column into two columns at the space
    df[['Sample', 'Count']] = df['Data'].str.split(' ', n=1, expand=True)

    # Extract relevant information
    sample_name = df.iloc[0, 0].split('_')[0]  # Extract 'IU112_S101'

    # Create a new DataFrame
    new_df = pd.DataFrame(index=[sample_name])

    # Iterate through rows and populate the new DataFrame
    for index, row in df.iterrows():
        column_name = row['Sample'].split('_')[-1].replace('.bam', '')  # Extract column name
        value = row['Count']
        new_df.loc[sample_name, column_name] = value

    # Remove the original column if needed
    df = df.drop(columns=['Data']) # Dropping the original column named 'Data'
    df = df.iloc[1:]

    # Extract "mouseBetter" using string manipulation
    df['Sample'] = df['Sample'].str.split('_').str[-1].str.split('.').str[0]

    # Store the DataFrame in the dictionary
    dfs[df_name] = df

# Access the DataFrames using their names
# For example, to access the DataFrame for "IU112_S101_read_counts.txt":
# df_IU112_S101 = dfs["IU112_S101"]

dfs['IU112_S101_read_counts.txt']

# List to hold individual DataFrames with 'Category' as index
dfs_with_index = []

for sample, df in dfs.items():
    # Set 'Sample' as index and rename the 'Count' column to the sample name
    df = df.set_index('Sample')
    df = df.rename(columns={'Count': sample})
    dfs_with_index.append(df)

# Merge all dataframes on the 'Category' index
merged_df = pd.concat(dfs_with_index, axis=1)

merged_df = merged_df.rename(columns=lambda x: x.replace('_read_counts.txt', '') if '_read_counts.txt' in x else x)

df=merged_df.T
df

# Assuming df is your DataFrame
df = df.drop(columns=['mouseLoss', 'humanLoss'])
df

# Assuming df is your DataFrame
df['humanBetter'] = df['humanBetter'].astype(int)
df['mouseBetter'] = df['mouseBetter'].astype(int)
df['humanOnly'] = df['humanOnly'].astype(int)
df['mouseOnly'] = df['mouseOnly'].astype(int)
df['human'] = df['humanBetter'] + df['humanOnly']
df['mouse'] = df['mouseBetter'] + df['mouseOnly']

# Optionally, you can remove the original columns if you no longer need them:
df = df.drop(columns=['humanBetter', 'humanOnly', 'mouseBetter', 'mouseOnly'])
df.info()

df

# Calculate total counts for each row
df['total'] = df['human'] + df['mouse']

# Calculate percentages for each column
df['human_pct'] = round((df['human'] / df['total']) * 100,2)
df['mouse_pct'] = round((df['mouse'] / df['total']) * 100,2)

# Create the new dataframe with only percentages
percentage_df = df[['human_pct', 'mouse_pct']]

# Display the new dataframe
percentage_df

df = df[['human','mouse']]
df

import matplotlib.pyplot as plt
import pandas as pd

# Create the stacked bar plot
ax = df.plot(kind='bar', stacked=True, figsize=(10, 6))

# Set plot labels and title
plt.xlabel('Sample')
plt.ylabel('Number of Aligned reads')
plt.title('Human-Aligned and Mouse-Aligned Reads')

# Rotate x-axis labels for better readability
plt.xticks(rotation=45, ha='right')

# Display the plot
plt.show()

!jupyter nbconvert --to html /content/your_notebook.ipynb
