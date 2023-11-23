import os
import pandas as pd
import re


run_name = "t1_231123"#input("What is the run name? ")
grd_name = "175-235"#input("The grid cell?lat-long ")
path_csv = f"/home/bianca/bianca/CAETE-DVM-alloc-allom/outputs/{run_name}/gridcell{grd_name}/csv"
start_year = "1979" #input("what is the start year? ")
end_year = "2017" #input("what is the start year? + 1")
start_year = int(start_year)
end_year = int(end_year)

# Get a list of all files with .csv extension in the directory(there is one csv for each alive PLS)
list_files = [file for file in os.listdir(path_csv) if file.endswith(".csv")]
   

# # # Extract the final group of numbers from the file names(get the PLS id)
s = [re.search(r"_([0-9]+)\.csv", file).group(1) for file in list_files if re.search(r"_([0-9]+)\.csv", file)]

# Create an empty DataFrame to store the merged and sorted data
final_merged_df = pd.DataFrame()

# Extract the final group of numbers from the file names and group the files (group each PLS considering all spins)
for file in list_files:
    match = re.search(r"EV_([0-9]+)\.csv", file)
    if match:
        group_number = match.group(1)
        file_path = os.path.join(path_csv, file)
        df = pd.read_csv(file_path)
        # Add a column to store the group number for reference
        df['GroupNumber'] = group_number
        # Merge the current group into the final DataFrame
        final_merged_df = pd.concat([final_merged_df, df], ignore_index=True)

# Create a DataFrame with all combinations of years and group numbers
all_years = range(start_year, end_year)

#Select the PLS ID for the alives
all_group_numbers = final_merged_df['GroupNumber'].unique()



all_combinations = pd.DataFrame([(year, group_number) for year in all_years for group_number in all_group_numbers],
                                 columns=['YEAR', 'GroupNumber'])



# Merge all_combinations with final_merged_df to fill gaps
final_merged_df = pd.merge(all_combinations, final_merged_df, on=['YEAR', 'GroupNumber'], how='left')

# Fill NaN values in 'PID' with 'GroupNumber' and convert back to int
final_merged_df['PID'] = final_merged_df['PID'].fillna(final_merged_df['GroupNumber']).astype(int)

final_merged_df['OC'] = final_merged_df['OC'].fillna(0.0)  # Fill NaN values in 'OC' with 0

# Sort the final DataFrame by the "YEAR" column in ascending order
final_merged_df.sort_values(by=["GroupNumber", "YEAR"], inplace=True)

# Drop the temporary group number column
final_merged_df.drop(columns=['GroupNumber'], inplace=True)

# Save the final merged and sorted DataFra  prinme to a new CSV file
final_file_path = os.path.join(path_csv, "final_merged_sorted_data.csv")
final_merged_df.to_csv(final_file_path, index=False)

#PLS id, year and occupation
PLS_ocp_year = final_merged_df

# Get a list of all files with .csv extension in the directory
list_files = [file for file in os.listdir(path_csv) if file.startswith("sorted_merged_") and file.endswith(".csv")]

# Delete temporary files
for file in list_files:
    file_path = os.path.join(path_csv, file)
    os.remove(file_path)

print("Your file has been created! Find it in:", path_csv)

#Now get the trait values from attrs (without considering the occupation)
# Read file with all pls traits
pls_traits = pd.read_csv("/home/bianca/bianca/CAETE-DVM-alloc-allom/outputs/pls_attrs-3000.csv")

# Get the PIDs, that is, the alive PLSs
pids_to_select = final_merged_df['PID'].unique()

# Filtrar pls_traits com base nos PIDs em pids_to_select
selected_PLS_traits = pls_traits[pls_traits['PLS_id'].isin(pids_to_select)]

# Realiza a agregação
ocp_traits = pd.merge(PLS_ocp_year, selected_PLS_traits[['PLS_id', 'sla_random']], left_on='PID', right_on='PLS_id', how='left')
# Remove a coluna 'PLS_id' da nova tabela
ocp_traits = ocp_traits.drop('PLS_id', axis=1)


#Calculates the value for a trait multiplying it by the PLS occupation
ocp_traits['sla_ocp'] = ocp_traits['OC']*ocp_traits['sla_random']

# Save the final merged and sorted DataFra  prinme to a new CSV file
final_merged_path = os.path.join(path_csv, "PLS_alive_traits.csv")
ocp_traits.to_csv(final_merged_path, index=False)






# # # # Select only data from living pls (those in the list of files)
# # new_pls_traits = pls_traits[pls_traits['PLS_id'].isin(pls_id['pls_id'])]
    

# # Merge the all_combinations DataFrame with the final DataFrame to fill gaps
# final_merged_df = pd.merge(all_combinations, final_merged_df, on=['YEAR', 'GroupNumber'], how='left')

# # Fill NaN values in columns other than 'YEAR' and 'GroupNumber' with 0
# final_merged_df.fillna(0, inplace=True)

# # Sort the final DataFrame by the "YEAR" column in ascending order
# final_merged_df.sort_values(by=["GroupNumber", "YEAR"], inplace=True)

# # Drop the temporary group number column
# final_merged_df.drop(columns=['GroupNumber'], inplace=True)

# # Save the final merged and sorted DataFrame to a new CSV file
# final_file_path = os.path.join(path_csv, "final_merged_sorted_data.csv")
# final_merged_df.to_csv(final_file_path, index=False)

# # # Sort the final DataFrame by the "YEAR" column in ascending order
# # final_merged_df.sort_values(by=["GroupNumber", "YEAR"], inplace=True)

# # # Drop the temporary group number column
# # final_merged_df.drop(columns=['GroupNumber'], inplace=True)

# # # Save the final merged and sorted DataFrame to a new CSV file
# # final_file_path = os.path.join(path_csv, "final_merged_sorted_data.csv")
# # final_merged_df.to_csv(final_file_path, index=False)



# # # Print the file groups
# # # for group_number, files in file_groups.items():
# #     # print(f"Files with group number {group_number}: {files}")






# # # # Create a DataFrame with the extracted numeric part
# # pls_id = pd.DataFrame({"pls_id": pd.to_numeric(s)})

# # # # Check the structure of pls_id to make sure it is as expected
# # print(pls_id.info())

# # # # Read file with all pls traits
# # pls_traits = pd.read_csv("/home/bianca/bianca/CAETE-DVM-alloc-allom/outputs/pls_attrs-3000.csv")

# # # # Select only data from living pls (those in the list of files)
# # new_pls_traits = pls_traits[pls_traits['PLS_id'].isin(pls_id['pls_id'])]

