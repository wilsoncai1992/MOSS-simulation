# Simulation of Cai W, van der Laan MJ (2019+). *One-step TMLE for time-to-event outcomes.*


## Pre-requisites

```
# MOSS package version >= 1.1.2
devtools::install_github("wilsoncai1992/MOSS")
install.packages("survtmle")
```


## Instructions

To reproduce the simulation section

```R
cd ./code_simulation/
# run simulation
R CMD BATCH ./run_simulation.R
mkdir ./output/
mv ./df_metric.rda ./output/df_metric.rda
# create plots
R CMD BATCH ./plot_result.R
# the plots will be saved in `./code_simulation/output/` after the script
```

To reproduce the real data analysis section

```R
cd ./code_mgus2/
R CMD BATCH ./run_simulation.R
# the plots will be saved in `./code_mgus2/` after the script
```

