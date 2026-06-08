# THAMES
This repository contains the code files for the THAMES project. The purpose of this project was to develop an accurate custom simulator for LLCs made of STTRAM in C++, combining accurate spin physics with high granularity compact RC modelling. The files in this repository can be used to obtain thermal modelling estimates for any STTRAM cache as long as a gem5 or similar CPU simulator byte-level output trace is available for the program you want to test.

# Files
- Read_sim_stats.h: reads gem5 byte-level trace data for a workload.
- STTRAM_Power_model.h: Contains bit and 4kiB-level power consumption calculations for MTJs of specific parameters (40nm MTJ, 22nm CMOS logic); change constants to suit the MTJ and circuitry you are modelling.
- STTRAM_Thermal_Model.cpp: Updates thermals for a 8MiB L2 cache using RC modelling. Only edit if you want to model a cache of a different size or want to change some of the material parameters being used (like the cooling rate of the heatsink or the size/material of the heatspreader/TIM).
- main.cpp: Executes the simulation. Only edit the line where the file path for reading the gem5 trace is defined do not touch anything else.
