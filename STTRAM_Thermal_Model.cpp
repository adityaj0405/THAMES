#include <vector>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <numeric>
#include <filesystem>
#include <array>
#include <bit>
#include "Read_sim_stats.h"
#include "STTRAM_Power_model.h"
#include "STTRAM_thermal_model.h"
#include <fstream>
namespace fs = std::filesystem;



// --- CONSTRUCTOR ---
ThermalModel::ThermalModel() : 
    R_vertical_bank(setResistance(k_si, t_si, A_bank)),
    R_lateral_bank(setResistance(k_si, t_si_lateral, A_bank_cross_section)),
    C_thermal_bank(setCapacitance(c_si, t_si, A_bank, scaling_factor_si)),
    R_vertical_beol(setResistance(k_beol, t_beol, A_bank)),
    C_thermal_beol(setCapacitance(c_beol, t_beol, A_bank, 1.0)),
    R_vertical_spreader(setResistance(k_cu, t_cu, A_spreader) + R_vertical_TIM),
    C_thermal_spreader(setCapacitance(c_cu, t_cu, A_spreader, scaling_factor_cu)) // Calculate R and C values for the banks and the spreader.
    {
    setPHoldAndInitialiseTemps(); //Calculate the power needed to hold the banks at 70C when the cache is idle
}


// --- HELPER FUNCTIONS ---
void ThermalModel::loadSimData(std::string path){
    sim_write_data = parseSimulatorCSV(path); //Load byte-level write traces from gem5's output
    current_workload_name = fs::path(path).stem().string(); //Get the name of the workload from the PARSEC suite
}


int ThermalModel::map2Dto1D(int i, int j){
    return (i * vert_dims) + j; // Convert a coordinate in the 2D grid of banks to a single value corresponding to a row in the data
}


double ThermalModel::setResistance(double k, double t, double A) {
    return t /(k * A);
}


double ThermalModel::setCapacitance(double c, double t, double A, double scaling_factor) {
    return c * t * A * scaling_factor;
}


void ThermalModel::setPHoldAndInitialiseTemps(){
    double delta_T = T_ss - T_ambient; 
    
    // Calculate the equivalent resistance of the entire package
    double R_total_package = (R_vertical_bank / N_banks) + R_vertical_spreader + R_vertical_sink;
    
    // Total Watts required for the whole chip to hold 70C
    double P_total_package = delta_T / R_total_package;
    
    // Divide into resistance per banl
    P_hold = P_total_package / N_banks;


    for (int i = 0; i < hor_dims; i++) {
        for (int j = 0; j < vert_dims; j++) {
            T_feols_current[i][j] = T_ss;
            T_feols_next[i][j] = T_ss;

            T_beols_current[i][j] = T_ss;
            T_beols_next[i][j] = T_ss;
            

            P_feols[i][j] = 0.0;
            P_beols[i][j] = 0.0;
        }
    }


    // Initialise the package nodes based on 70 degree target for silicon and ambient temps
    //calculate their exact resting temperatures based on the heat flowing through them
    T_sink_current = T_ambient + (P_total_package * R_vertical_sink);
    T_sink_next = T_sink_current;
    
    T_spreader_current = T_sink_current + (P_total_package * R_vertical_spreader);
    T_spreader_next = T_spreader_current;
}


void ThermalModel::resetPower(){
    for (int i = 0; i < hor_dims; i++){
        for (int j = 0; j < vert_dims; j++){
            P_feols[i][j] = 0.0;
            P_beols[i][j] = 0.0;
        }
    }
}



// Calculate power dissipated in each bank at the current timestep
void ThermalModel::getPower(int i, int j, int timestep) {
    double T_current_bank = T_feols_current[i][j]; // Get temperature of (i, j)-th bank
    int row_idx = map2Dto1D(i, j); // Convert 2D bank location to 1D location in csv file
    int N_writes = sim_write_data[row_idx][timestep]; //Parse the number of writes at the timestep
    if (N_writes == 0){
        return;
    }
    std::pair<double, double> E_bank = calculateTotalBankEnergy(T_current_bank, N_writes, beta); // Calculate total energy dissipated in bank
    double P_feol = E_bank.second/dt_csv + P_leakage; // Divide by timestep to get averaged feol power from peripherals and ewt, add leakage on if bank not gated
    double P_beol = E_bank.first/dt_csv;

    if (i >= 0 && i < hor_dims && j >= 0 && j < vert_dims) {
        P_feols[i][j] += P_feol; //Check to make sure i and j are within the size of the banks and push the power consumption to the power arrays
        P_beols[i][j] += P_beol;
    }
}


void ThermalModel::logTimestepData(std::ofstream& log_file, int timestep) {
    for (int i = 0; i < hor_dims; i++) {
        for (int j = 0; j < vert_dims; j++) {
            // Convert to Celsius for logging
            double temp_beol_C = T_beols_current[i][j] - 273.15; 
            double temp_feol_C = T_feols_current[i][j] - 273.15; 
            
            log_file << timestep << "," 
                     << i << "," 
                     << j << "," 
                     << P_beols[i][j] << "," 
                     << P_feols[i][j] << "," 
                     << temp_beol_C << ","
                     << temp_feol_C << "\n";
        }
    }
}


// --- Lateral Heat Flow---
double ThermalModel::calculateLateralHeat(int i, int j) {
    double heat_loss = 0.0;
    double current_temp = T_feols_current[i][j];

    // Check left (i - 1)
    if (i > 0) heat_loss += (current_temp - T_feols_current[i-1][j]) / R_lateral_bank;
    // Check right (i + 1)
    if (i < 63) heat_loss += (current_temp - T_feols_current[i+1][j]) /R_lateral_bank;
    // Check top (j - 1)
    if (j > 0) heat_loss += (current_temp - T_feols_current[i][j-1]) /R_lateral_bank;
    // Check bottom (j + 1)
    if (j < 31) heat_loss += (current_temp - T_feols_current[i][j+1]) /R_lateral_bank;

    return heat_loss;
}

// --- RC equations to update temps ---
void ThermalModel::thermalStep() {
    double total_heat_to_spreader = 0.0;

    // EQUATION 1: Update the 2048 Subarrays
    for (int i = 0; i < hor_dims; i++) {
        for (int j = 0; j < vert_dims; j++) {

            // Vertical Heat Flow: BEOL -> Silicon
            // This represents heat moving through the 0.5um BEOL layer containing MTJs and the peripheral circuitry
            double heat_beol_to_si = (T_beols_current[i][j] - T_feols_current[i][j]) / R_vertical_beol;
            
            // Update BEOL Temperature (The MTJ layer)
            // Note: Dynamic power injected into the BEOL layer
            T_beols_next[i][j] = T_beols_current[i][j] + 
                (dt_sim / C_thermal_beol) * (P_beols[i][j]- heat_beol_to_si);
            
            double lateral_heat = calculateLateralHeat(i, j);
            double vertical_heat = (T_feols_current[i][j] - T_spreader_current) /R_vertical_bank;
            
            // Add to the total going up to the spreader for Eq 2
            total_heat_to_spreader += vertical_heat;

            // Calculate the new temperature and store it in the next array
            T_feols_next[i][j] = T_feols_current[i][j] + 
                (dt_sim / C_thermal_bank) * (P_hold + P_feols[i][j] - lateral_heat - vertical_heat + heat_beol_to_si);
            
            //T_beols_next[i][j] = T_feols_next[i][j];
        }
    }

    // EQUATION 3: Update the Spreader
    double heat_to_sink = (T_spreader_current - T_sink_current) / R_vertical_spreader;
    T_spreader_next = T_spreader_current + 
        (dt_sim / C_thermal_spreader) * (total_heat_to_spreader - heat_to_sink);

    // EQUATION 4: Update the Sink
    double heat_to_ambient = (T_sink_current - T_ambient) / R_vertical_sink;
    T_sink_next = T_sink_current + 
        (dt_sim / C_thermal_sink) * (heat_to_sink - heat_to_ambient);


    /* Overwrite current arrays with next arrays only after the next arrays 
       have been fully calculated using current information */
    

    for (int i = 0; i < hor_dims; i++) {
        for (int j = 0; j < vert_dims; j++) {
            T_beols_current[i][j] = T_beols_next[i][j];
            T_feols_current[i][j] = T_feols_next[i][j];
        }
    }

    T_spreader_current = T_spreader_next;
    T_sink_current = T_sink_next;
}



void ThermalModel::runFullSimulation(){
    int N_timesteps = sim_write_data[0].size();

    // DEBUG: Find the exact running location of the program
    std::cout << "DEBUG: Current working directory is: " << fs::current_path().string() << "\n";
    // Define the results folder wherever is most convenient for you
    std::string results_folder = "C:/Users/Aditya/Documents/Uni/Individual Project/Results";
    std::error_code ec; //This will catch OS errors
    
    // Check if folder exists, and catch any OS errors
    if (!fs::exists(results_folder, ec)) {
        // Try to create the directory, catch OS errors if it fails
        if (!fs::create_directory(results_folder, ec)) {
            std::cerr << "CRITICAL ERROR: OS refused to create 'Results' folder!\n";
            std::cerr << "Reason: " << ec.message() << "\n";
            return; // Abort the simulation if the results folder cant be created
        }
    }

    // Build the clean log filename
    std::string out_filename = current_workload_name + "_thermal_output.csv";
    fs::path full_path = fs::path(results_folder) / out_filename;
    
   
    std::cout << "DEBUG: Attempting to save to: " << fs::absolute(full_path).string() << "\n";

    // Open the log file
    std::ofstream log_file(full_path);
    
    if (log_file.is_open()) {
        std::cout << "SUCCESS: Log file created and opened!\n";
        log_file << "Timestep,Bank_X,Bank_Y,Power_BEOL, Power_FEOL, Temp_BEOL, Temp_FEOL\n"; 
    } else {
        std::cerr << "CRITICAL ERROR: Failed to open file at " << full_path.string() << "\n";
        std::cerr << "Workload name parsed as: '" << current_workload_name << "'\n"; 
        return; 
    }

    for (int timestep = 0; timestep < N_timesteps; timestep++){
        for (int i = 0; i < 64; i++) {
            for (int j = 0; j < 32; j++) {
                getPower(i, j, timestep); //Increment timestep and compute the power dissipation per bank 
            }
        }

        for (int subtick = 0; subtick < subticks; subtick++){
            thermalStep(); //thermal step for however many subticks for numerical stability rather than trying to update 33us worth of temp updates at once
        }

        if (log_file.is_open()) {
            logTimestepData(log_file, timestep);
        }

        resetPower(); //Reset banks to idle/static leakage power

        if (timestep  % 5 == 0) {
            printMaxTemp();// For debugging
        }
    }
}



// --- DEBUGGING ---
void ThermalModel::printMaxTemp() {
    double max_temp_beol = 0.0;
    double max_temp_feol = 0.0;
    int max_i = 0;
    int max_j = 0;

    // Scan all 2048 banks for the hottest temperature
    for (int i = 0; i < hor_dims; i++) {
        for (int j = 0; j < vert_dims; j++) {
            if (T_beols_current[i][j] > max_temp_beol) {
                max_temp_beol = T_beols_current[i][j];
                max_temp_feol = T_feols_current[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    std::cout << "Center Bank [32][16] Temp: " << T_beols_current[32][16] - 273.15 << " C\n";
    std::cout << "Hottest Bank BEOL Temp [" << max_i << "][" << max_j << "] Temp: " << max_temp_beol - 273.15 << " C\n";
    std::cout << "Hottest Bank FEOL Temp [" << max_i << "][" << max_j << "] Temp: " << max_temp_feol - 273.15 << " C\n";
    std::cout << "Spreader Temp:             " << T_spreader_current - 273.15 << " C\n";
    std::cout << "Sink Temp:                 " << T_sink_current - 273.15 << " C\n";
    std::cout << "-----------------------------------\n";
}
