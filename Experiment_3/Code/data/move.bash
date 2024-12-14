for file in PHYS506_Experiment_3_Data*.csv; do
  # Extract the part between parentheses
  new_name=$(echo "$file" | sed -E 's/PHYS506_Experiment_3_Data\((.*)\)\.csv/\1.csv/')
  mv "$file" "$new_name"
done
