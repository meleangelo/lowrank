---
output:
  pdf_document: default
  html_document: default
---
# Low-rank simulations

This repo contains simulations to check the low-rank approximation of conditional probabilities in the dynamic network model.

The repo contains several folders:

- `code`: contains all code to run the simulations
- `output`: simulation output is saved and stored here
- `notes`: latex notes (includes graphs and tables)


## `code`
Each simulation needs 3 files. 

- `setup_functions.R` contains all the functions
- `simulation_design#.R` contains the design (see table below for the complete list)
- `simulate_loop.R` runs the simulation, saves results, and outputs a .tex table.

All the simulations are contained in the file `simulation.R`. 

### simulation designs

Design | $n$          | $T$ | $K$ | $d$ | $\gamma$ | $\nu$
:-:    | :-:          | :-: | :-: | :-: | :-:      | :-: 
1      | 2000-10000   | 5   | 2   | 1   | 0.1      | (0.8, -1.5)
2      | 2000-10000   | 5   | 2   | 1   | 0.5      | (0.8, -1.5)
3      | 2000-10000   | 5   | 2   | 1   | 1.0      | (0.8, -1.5)
4      | 2000-10000   | 5   | 2   | 1   | 1.5      | (0.8, -1.5)
5      | 2000-10000   | 5   | 4   | 1   | 0.1      | (0.8, -1.5)
6      | 2000-10000   | 5   | 4   | 1   | 0.5      | (0.8, -1.5)
7      | 2000-10000   | 5   | 4   | 1   | 1.0      | (0.8, -1.5)
8      | 2000-10000   | 5   | 4   | 1   | 1.5      | (0.8, -1.5)


TO-DO
