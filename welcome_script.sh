#!/bin/bash

# Colors for pretty output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to clear screen
clear_screen() {
clear
}

# Function to check if python and required packages are installed
check_sims() {
echo -e "${BLUE}Checking simulation directories...${NC}"
python src/check_sims.py
echo -e "\nPress Enter to return to main menu..."
read
}

# Function to compute power spectrum
compute_power_spectrum() {
echo -e "${BLUE}=== Power Spectrum Computation ===${NC}"

# Get simulation path using select
echo -e "\nAvailable simulation directories:"
echo -e "${GREEN}Please select a directory:${NC}"

# Create array of directories
directories=($(find data -maxdepth 2 -type d | sort))

# Check if any directories were found
if [ ${#directories[@]} -eq 0 ]; then
    echo -e "${RED}No simulation directories found in data/!${NC}"
    echo -e "\nPress Enter to return to main menu..."
    read
    return
fi

select simulation_path in "${directories[@]}" "Return to main menu"; do
    if [ "$simulation_path" = "Return to main menu" ]; then
        return
    elif [ -n "$simulation_path" ]; then
        break
    else
        echo -e "${RED}Invalid selection. Please try again.${NC}"
    fi
done

# Check if path exists (redundant but safe)
if [ ! -d "$simulation_path" ]; then
    echo -e "${RED}Error: Directory '$simulation_path' does not exist!${NC}"
    echo -e "\nPress Enter to return to main menu..."
    read
    return
fi

# Get optional parameters with defaults
echo -e "\nEnter number of mesh cells [default: 256]:"
read nmesh
nmesh=${nmesh:-256}

echo -e "Enter box size in Mpc/h [default: 100.0]:"
read box_size
box_size=${box_size:-100.0}

echo -e "\n${GREEN}Computing power spectrum with:${NC}"
echo "Simulation path: $simulation_path"
echo "Nmesh: $nmesh"
echo "Box size: $box_size Mpc/h"

# Run the Python script
python src/main.py --input_path "$simulation_path" --nmesh "$nmesh" --box-size "$box_size"

echo -e "\nPress Enter to return to main menu..."
read
}

# Function to generate comparison plots
generate_comparison() {
echo -e "${BLUE}=== Simulation Comparison Plots ===${NC}"

# Get optional parameters with defaults
echo -e "\nEnter box size in Mpc/h [default: 100.0]:"
read box_size
box_size=${box_size:-100.0}

echo -e "Enter number of mesh cells [default: 256]:"
read nmesh
nmesh=${nmesh:-256}

echo -e "\n${GREEN}Generating comparison plots...${NC}"
python src/LR_HR_SR.py --box-size "$box_size" --nmesh "$nmesh"

echo -e "\nPress Enter to return to main menu..."
read
}

# Main menu loop
while true; do
clear_screen
echo -e "${BLUE}=== Power Spectrum Analysis Tool ===${NC}"
echo -e "\nPlease choose an option:"
echo "1. Run Environment Check"
echo "2. Compute Power Spectrum"
echo "3. Generate Comparison Plots"
echo "4. Exit"

echo -e "\nEnter your choice (1-4):"
read choice

case $choice in
    1)
        check_sims
        ;;
    2)
        compute_power_spectrum
        ;;
    3)
        generate_comparison
        ;;
    4)
        echo -e "\n${GREEN}Exiting program. Goodbye!${NC}"
        exit 0
        ;;
    *)
        echo -e "${RED}Invalid option. Please try again.${NC}"
        sleep 1
        ;;
esac
done