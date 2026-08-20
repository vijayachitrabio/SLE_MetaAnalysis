import pandas as pd

# Read the bivariate results
df = pd.read_csv("results/lava_crosstrait_bivariate.csv")

# Filter for tests involving SLE_Discovery
df_sle = df[(df['phen1'] == 'SLE_Discovery') | (df['phen2'] == 'SLE_Discovery')].copy()

# Ensure SLE is always phen1 for consistency
for idx, row in df_sle.iterrows():
    if row['phen2'] == 'SLE_Discovery':
        df_sle.at[idx, 'phen1'] = 'SLE_Discovery'
        df_sle.at[idx, 'phen2'] = row['phen1']

# Keep relevant columns and rename for publication
# LAVA columns: phen1, phen2, rho, rho.lower, rho.upper, r2, r2.lower, r2.upper, p, LOC
df_sle = df_sle[['LOC', 'phen2', 'rho', 'rho.lower', 'rho.upper', 'p']]
df_sle.rename(columns={
    'phen2': 'Secondary_Trait',
    'rho': 'Genetic_Correlation_Rho',
    'rho.lower': 'Rho_95CI_Lower',
    'rho.upper': 'Rho_95CI_Upper',
    'p': 'P_value'
}, inplace=True)

# Format the floats
df_sle['Genetic_Correlation_Rho'] = df_sle['Genetic_Correlation_Rho'].round(3)
df_sle['Rho_95CI_Lower'] = df_sle['Rho_95CI_Lower'].round(3)
df_sle['Rho_95CI_Upper'] = df_sle['Rho_95CI_Upper'].round(3)
df_sle['P_value'] = df_sle['P_value'].apply(lambda x: f"{x:.2e}")

# Read the loci information to get gene names
loci_df = pd.read_csv("results_extracted/final_lava_consolidated_loci.tsv", sep="\t")
loci_df = loci_df[['LAVA_LOC', 'Gene', 'Region', 'RSID']]
loci_df.rename(columns={'LAVA_LOC': 'LOC'}, inplace=True)

# Merge
final_df = pd.merge(df_sle, loci_df, on='LOC', how='left')

# Reorder columns
final_df = final_df[['LOC', 'Region', 'Gene', 'RSID', 'Secondary_Trait', 'Genetic_Correlation_Rho', 'Rho_95CI_Lower', 'Rho_95CI_Upper', 'P_value']]

# Sort by p-value
# final_df['P_value_float'] = final_df['P_value'].astype(float)
# final_df = final_df.sort_values('P_value_float').drop(columns=['P_value_float'])

# Save to TSV
assert final_df.duplicated(subset=['LOC', 'Secondary_Trait']).sum() == 0, "Duplicate LOC-Trait pairs found!"
final_df.to_csv("results/lava_crosstrait_results.tsv", sep="\t", index=False)
print("Saved final results to results/lava_crosstrait_results.tsv")
