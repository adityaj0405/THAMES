#include <iostream>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>
#include <array>
#include <bit>

constexpr double GAMMA_GYRO     = 1.76086e11;   // rad/(s*T)
constexpr double MU0            = 1.2566e-6;    // T*m/A

// --- 40nm p-MTJ Properties (Sources: Ikeda et al, 2010, Yakushiji et. al 2010) ---
constexpr double ALPHA          = 0.027;        // Gilbert Damping
constexpr double MU0_H_EFF      = 0.28;         // Tesla (Effective Anisotropy)
constexpr double IC0_AT_300K    = 49.0e-6;      // Amps (Critical Current @ 300K)
constexpr double DELTA_THERMAL  = 37.6;         // Thermal Stability Factor
constexpr double d_nm           = 40;   // MTJ Diameter (40nm)
constexpr double r_μm           = 0.02; //MTJ Radius in μm
const double A = M_PI * pow(r_μm, 2);
constexpr double RA             = 4.4;          // Resistance-Area Product Ωμm^2
constexpr double TMR            = 0.84;         // TMR Ratio (R_AP/R_A)

// Initial Angle derived from Delta: sqrt(1 / 2*Delta)
const double THETA_0 = std::sqrt(1.0 / (2.0 * DELTA_THERMAL)); 

// --- Thermal Degradation (Bloch Law Linearization) ---
constexpr double K_FIT          = 0.0026;        // 1/K 

// --- Circuit / Architecture (Source: NVSim / ITRS 22nm) ---
constexpr double V_DD           = 0.9;               // Volts
constexpr double R_ACCESS       = 4000.0;            // Ohms (Access Transistor)
const double R_P                = RA/A;              // Ohms (Parellel MTJ Resistance))
const double R_AP               = R_P * (1 + TMR);   // Ohms (Antiparallel MTJ Resistance)
const double R_MTJ_AVG          = (R_AP + R_P)/2;    // Ohms (Average of P and AP States)
constexpr double R_BL           = 500;               // Ohms (Bitline Resistance)
constexpr double E_peripheral   = 0.0802e-9;         // Zhou et. al peripheral overhead scaled to 22nm
constexpr double E_EWT          = 0.0181e-9;         // Zhou et. al EWT overhead scaled to 22nm
constexpr double CACHE_LINE_SIZE = 64.0;

// --- Helper Functions ---
const double I_WRITE            = (V_DD/(R_MTJ_AVG+R_ACCESS+R_BL));     // Amps (Write Current calculated from path resistance and supply voltage)

inline double calculateCriticalCurrent(double T){
    double IC0_AT_TEMP_T = IC0_AT_300K*(1-K_FIT*(T-300));
    return IC0_AT_TEMP_T; // Linear approximation of Bloch's law (see paper)
}

inline double calculateWriteTime(double T){
    double I_C0_AT_TEMP_T = calculateCriticalCurrent(T);
    double overdrive = I_WRITE/I_C0_AT_TEMP_T;
    double t_char = (log(M_PI/(2*THETA_0)))/(ALPHA*GAMMA_GYRO*MU0_H_EFF);
    double t_write = t_char/(overdrive-1);
    return t_write; // Solution to LLG for average switching time (see paper)
}

inline double calculateWriteEnergy(double T, double t_write) {
    double E_write = V_DD * I_WRITE * t_write;
    return E_write; 
}




inline std::pair<double, double> calculateTotalBankEnergy(double T_bank, int N_writes, double beta) {
    double N_accesses = std::ceil(N_writes / CACHE_LINE_SIZE); //Overhead energy in Zhou et. al is given per access, so scale no. of bytes written to number of cache lines 
    double E_total = 0.0;
    std::pair<double, double> Energies;
    double t_writes = calculateWriteTime(T_bank);
    double E_writes = calculateWriteEnergy(T_bank, t_writes);
    double total_write_energy = N_writes * 8 * beta * E_writes;
    double total_overhead = (E_EWT + E_peripheral) * N_accesses;
    Energies.first = total_write_energy;
    Energies.second = total_overhead;
    return Energies;
}

