import os
import pandas as pd
import re


run_name = "tresmil"#input("What is the run name? ")
grd_name = "175-235"#input("The grid cell?lat-long ")
path_csv = f"/home/bianca/bianca/CAETE-DVM-alloc-allom/outputs/{run_name}/gridcell{grd_name}/csv"

# Get a list of all files with .csv extension in the directory
list_files = [file for file in os.listdir(path_csv) if file.endswith(".csv")]
   

# # # Extract the final group of numbers from the file names
s = [re.search(r"_([0-9]+)\.csv", file).group(1) for file in list_files if re.search(r"_([0-9]+)\.csv", file)]

# Create an empty DataFrame to store the merged and sorted data
final_merged_df = pd.DataFrame()

# Extract the final group of numbers from the file names and group the files
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

# Sort the final DataFrame by the "YEAR" column in ascending order
final_merged_df.sort_values(by=["GroupNumber", "YEAR"], inplace=True)

# Drop the temporary group number column
final_merged_df.drop(columns=['GroupNumber'], inplace=True)

# Save the final merged and sorted DataFrame to a new CSV file
final_file_path = os.path.join(path_csv, "final_merged_sorted_data.csv")
final_merged_df.to_csv(final_file_path, index=False)



# Print the file groups
# for group_number, files in file_groups.items():
    # print(f"Files with group number {group_number}: {files}")






# # # Create a DataFrame with the extracted numeric part
# pls_id = pd.DataFrame({"pls_id": pd.to_numeric(s)})

# # # Check the structure of pls_id to make sure it is as expected
# print(pls_id.info())

# # # Read file with all pls traits
# pls_traits = pd.read_csv("/home/bianca/bianca/CAETE-DVM-alloc-allom/outputs/pls_attrs-3000.csv")

# # # Select only data from living pls (those in the list of files)
# new_pls_traits = pls_traits[pls_traits['PLS_id'].isin(pls_id['pls_id'])]

