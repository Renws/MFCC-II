#!/bin/bash

start=$EPOCHREALTIME

# ========================= Configuration Section =========================
# Set the target directory prefix (e.g., mol00, mol01, mol02)
# You can modify this value manually, or pass it as a command line argument
TARGET_DIR_PREFIX="mol00"  # Default value: mol00
# =========================================================================

# Function: Show usage instructions
show_usage() {
    echo "Usage: $0 [directory_prefix]"
    echo "Example: $0 mol00  # Calculate energy for mol00 directory"
    echo "         $0 mol01  # Calculate energy for mol01 directory"
    exit 1
}

# Check command line arguments (optional: override default directory prefix)
if [ $# -eq 1 ]; then
    TARGET_DIR_PREFIX="$1"
elif [ $# -gt 1 ]; then
    echo "Error: Too many arguments!" >&2
    show_usage
fi

# Define subdirectories (fixed: tmpfile and dimertmpfile)
SUB_DIR1="tmpfile"          # Fixed subdir for 1B energy
SUB_DIR2="dimertmpfile"     # Fixed subdir for 2B energy

# Construct full commands with dynamic directory prefix
CMD1="python get_energy_1B.py --dir ${TARGET_DIR_PREFIX}/${SUB_DIR1}/"
CMD2="python get_energy_2B.py --dir ${TARGET_DIR_PREFIX}/${SUB_DIR2}/"

# Function: Extract the last column number from command output
extract_last_number() {
    local cmd="$1"
    # Execute the command and extract the last column number (only take the last line for multi-line output)
    local number=$(eval "$cmd" | tail -n 1 | awk '{print $NF}')
    # Validate if the extracted result is a valid number
    if ! [[ "$number" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
        echo "Error: Invalid number extracted from command '$cmd', extracted value: $number" >&2
        exit 1
    fi
    echo "$number"
}

# Extract the last column numbers from the output of the two programs
#echo "=== Processing directory: ${TARGET_DIR_PREFIX} ==="
#echo "Executing the first energy calculation program..."
num1=$(extract_last_number "$CMD1")
#echo "Number from first program (1B_energy): $num1"

#echo "Executing the second energy calculation program..."
num2=$(extract_last_number "$CMD2")
#echo "Number from second program (2B_energy): $num2"

# Perform high-precision floating-point calculation with bc, keep 6 decimal places
# Set scale=10 to calculate with 10 decimal places first, then format to 6 decimal places with printf to avoid precision loss
sum=$(echo "scale=10; $num1 + $num2" | bc)
formatted_sum=$(printf "%.6f" "$sum")

# Output the final result
echo -e "\n=== Calculation Result for ${TARGET_DIR_PREFIX} ==="
echo "1B_energy + 2B_energy = $formatted_sum"
echo "(energy_tot: $sum Hartree)"
end=$EPOCHREALTIME
awk -v s="$start" -v e="$end" 'BEGIN{printf "run_time：%.3fs\n", (e-s)}'
