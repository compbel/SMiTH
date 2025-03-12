# Loop from 1 to 500
for i in {1..500}
do
    echo "Processing iteration: $i"

    # Define paths for required input files
    patient_tree="pruned_trees_log50/FAVITES_output_logSI_contemp_T200_N100_E1_$i/patientTree.txt"
    patient_labeling="pruned_trees_log50/FAVITES_output_logSI_contemp_T200_N100_E1_$i/patientLabeling.txt"
    coloring="pruned_trees_log50/FAVITES_output_logSI_contemp_T200_N100_E1_$i/coloring.txt"
    output_dir="machina_results_with_CorrectRoot/FAVITES_output_logSI_contemp_T200_N100_E1_$i"
    root_file="Favites_log/FAVITES_output_logSI_contemp_T200_N100_E1_$i/error_free_files/transmission_network.txt"
    # Check if all required input files exist before proceeding
    if [[ ! -f "$patient_tree" || ! -f "$patient_labeling" || ! -f "$coloring" ]]; then
        echo "Skipping iteration $i: Missing required files."
        continue  # Skip this iteration
    fi

    # Create the output directory only if the files exist
    mkdir -p "$output_dir"

    # Extract the root number from transmission network.txt
    p_value=$(head -n 1 "$root_file" | awk '{print $2}')

    # Check if p_value is extracted correctly
    if [[ -z "$p_value" ]]; then
        echo "Error: Could not extract a valid -p value for iteration $i"
        continue  # Skip this iteration if extraction fails
    fi

    echo "Using -p $p_value for iteration $i"

    # Run the pmh command with the extracted value
    pmh -p "$p_value" -c "$coloring" "$patient_tree" "$patient_labeling" \
        -o "$output_dir" > "$output_dir/result.txt"

    echo "Finished iteration: $i"
done

echo "All iterations completed!"
