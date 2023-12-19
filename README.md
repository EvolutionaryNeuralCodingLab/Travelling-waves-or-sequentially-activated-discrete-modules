# Travelling-waves-or-sequentially-activated-discrete-modules

Repository for analyzing simulated or recorded electrophisiological data.

## General description
This code was used to analyze data for the paper Travelling-waves-or-sequentially-activated-discrete-modules by Orsher et al., 2023.

There are two main sections - **SimulationAnalysis** which contains the code used to analyze simualted 
data created by BrainStorm, and **WaveAnalysis** which contains the code used to
analyze the electrophysiological recordings (and some of the Simulated data).

To understand how the code was used, see the scripts in **analysis_scripts**, 
which show how the datas for the main figures were created.

Running analysis pipelines requires raw data which could be provided, upon a reasonable request, by [the Evolutionary Neural Coding Lab](https://www.evolutionaryneuralcodinglab.sites.tau.ac.il/).

## Dependencies
Data analysis of recorded electrophysiological uses [time-series-viewer](https://github.com/EvolutionaryNeuralCodingLab/time-series-viewer) code for data processing.
Please clone it and have add it to path.
