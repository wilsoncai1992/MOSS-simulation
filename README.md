# Simulation of Cai W, van der Laan MJ (2019+). One-step TMLE for time-to-event
outcomes.



## Instructions

To reproduce the simulation section

```R
R CMD BATCH ./R-paper/run_simulation.R
mkdir ./R-paper/output/
mv ./R-paper/df_metric.rda ./R-paper/output/df_metric.rda
# create plots
R CMD BATCH ./R-paper/plot_result.R
# the plots will be saved in `./R-paper/output/` after the script
```

To reproduce the real data analysis section

```
R CMD BATCH ./code_mgus2/run_simulation.R
# the plots will be saved in `./code_mgus2/` after the script
```

