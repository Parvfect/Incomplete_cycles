
import pandas as pd

results_filepath = r"C:\Users\Parv\Doc\RA\Projects\incomplete_cycles\data\clustering_post_synthesis_isolated\run_info_file.txt"


with open(results_filepath, 'r') as f:
    lines = f.readlines()

capping_dicts = []
no_capping_dicts = []

for line in lines:
    try:
        result_dict = eval(line)
        if result_dict['type'] == 'capping':
            capping_dicts.append(result_dict)
        else:
            no_capping_dicts.append(result_dict)
    except:
        continue

capping_df = pd.DataFrame(capping_dicts)
no_capping_df = pd.DataFrame(no_capping_dicts)

# Mean over the same coupling rates
capping_df_split = capping_df[['coupling_rate', 'best_recovery_clustering', 'mean_length', 'max_length', 'std_length', 'mean_deletions']]
no_capping_df_split = no_capping_df[['coupling_rate', 'best_recovery_clustering', 'mean_length', 'max_length', 'std_length', 'mean_deletions']]

grouped_capping = capping_df_split.groupby(['coupling_rate']).mean()
grouped_no_capping = no_capping_df_split.groupby(['coupling_rate']).mean()