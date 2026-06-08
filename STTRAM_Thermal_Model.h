#include <stdio.h>
#include <vector>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <numeric>
#include <array>
#include <bit>
#include <fstream>



class ThermalModel {
private:
    // State Variables (Temperature and Power)
    double T_feols_current[64][32];
    double T_feols_next[64][32];

    double T_beols_current[64][32];
    double T_beols_next[64][32];

    double P_feols[64][32];
    double P_beols[64][32];

    double T_spreader_current;
    double T_spreader_next;

    double T_sink_current;
    double T_sink_next;

    // Modelling Parameters
    static constexpr double A_bank = 1267.28e-12;
    static constexpr double w_bank = 35.6e-6;
    static constexpr int N_banks = 2048;
    static constexpr int hor_dims = 64;
    static constexpr int vert_dims = 32;
    static constexpr double T_ss = 343.15; //Steady-state target temperature
    static constexpr double T_ambient = 320.15; //Ambient temperature of case
    static constexpr double dt_csv = 0.000033; // Simulation timestep
    static constexpr double dt_sim = 1e-10;
    static constexpr int subticks = 330000;
    
    static constexpr double t_beol = 0.5e-6;
    static constexpr double k_beol = 1.4;
    static constexpr double c_beol = 1.96e6;



    static constexpr double t_si = 150e-6;
    static constexpr double t_si_lateral = w_bank; //Horizontal thickness of silicone heat has to travel through for centre of bank to centre of adjacent banks
    static constexpr double h_si_lateral = 51.38e-6;
    static constexpr double A_bank_cross_section = t_si_lateral * h_si_lateral;
    static constexpr double k_si = 120;
    static constexpr double c_si = 1.64e6;
    static constexpr double scaling_factor_si = 0.47;

    static constexpr double t_cu = 1e-3;
    static constexpr double k_cu = 390;
    static constexpr double c_cu = 3.47e6;
    static constexpr double scaling_factor_cu = 0.62;
    static constexpr double A_spreader = 9 * N_banks * A_bank;

    static constexpr double P_leakage = 64.7e-6;

    static constexpr double beta = 1.0;


    
    const double C_thermal_bank;
    const double R_vertical_bank;
    const double R_lateral_bank;

    const double C_thermal_beol;
    const double R_vertical_beol;

    constexpr static double R_vertical_TIM = 0.25;

    const double C_thermal_spreader;
    const double R_vertical_spreader;

    static constexpr double C_thermal_sink = 140.4;
    static constexpr double R_vertical_sink = 0.1;//Thermal Properties of Bank, Spreader and Sink materials

    double P_hold;


    std::vector<std::vector<double>> sim_write_data;
    std::string current_workload_name;


    






public:
    ThermalModel();

    int map2Dto1D(int i, int j);
    void loadSimData(std::string path);
    double setResistance(double k, double t, double A);
    double setCapacitance(double c, double t, double A, double scaling_factor);
    void setPHoldAndInitialiseTemps();


    void getPower(int i, int j, int timestep);
    void resetPower();
    double calculateLateralHeat(int i, int j);
    void thermalStep();
    void runFullSimulation();

 
    // For Testing Purposes
    void printMaxTemp();
    void logTimestepData(std::ofstream& log_file, int timestep);
};