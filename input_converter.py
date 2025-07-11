import pandas as pd

# Load the original CSV
df = pd.read_csv("sample_details_remove_incomplete_asp_glu.csv")

# Create the 'Protein' column
df['Protein'] = df['protein'].astype(str) + "_" + df['chain'].astype(str) + "_" + df['drug_id'].astype(str) + "_" + df['drug_name'].astype(str)

# Create a new DataFrame with only the 'Protein' column
output_df = df[['Protein']]

# Save to new CSV
output_df.to_csv("output.csv", index=False)

print("New CSV with 'Protein' column saved as 'output.csv'")
