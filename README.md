# 1D-Dyad-Model-RyR-IP3R
## Single dyad model configuration
### LTCC-triggered spark simulations
1. Perform LTCC-triggered spark simulations
    * Run: one_dyad_Mult_clusters_mult_runs_save_RyR_IP3R_6State.m (with IP3R model) or;
    * Run: one_dyad_Mult_clusters_mult_runs_save_RyR_ConstJSR.m (constant JSR leak) or;
    * Run: one_dyad_Mult_clusters_mult_runs_save_RyR_StochJSR.m (randomly occurring JSR leak) or;
    * Run: one_dyad_Mult_clusters_mult_runs_save_rand_RyR_IP3R.m (random cluster arrangement in dyad)
    * Run: one_dyad_Mult_clusters_mult_runs_save_RyR_Sheep.m (Sheep RyR, varying SERCA activity, no IP3R)
    * With dependencies:
        * RxnDiffusion_parameters.m
        * IP3R_parameters.m
        * EulerMethod.m
        * RxnDiffusion_Fluxes.m
        * Update_IP3R_vars.m
        * parsave.m
2. Result analysis & figure generation
    * Run: plot_Triggered_Spark_Cai_profile.m (temporal evolution of: dyad Ca<sup>2+</sup>, JSR Ca<sup>2+</sup>, n(RyR open), n(IP3R open)) 
    * With dependencies:
        * data_shade.m
    * Run: extract_Triggered_spark_data.m (identifies and extracts spark statistics into one file)
    * Run: plot_Triggered_Spark_Stats.m (swarm plots of spark amp and FDHM, statistical analysis)
    * With dependencies:
        * plotSpread.m
        * myErrorbar.m
        * isEven.m
        * repeatEntries.m

### Spontaneous spark simulations
1. Perform spontaneous spark simulations
    * Run: Coupl_rep_uncoupl_Mult_clusters_mult_runs_save_1Dyad_NoStim.m (with IP3R model) or;
    * Run: Coupl_rep_uncoupl_Mult_clusters_mult_runs_save_1Dyad_ConstJSR.m (constant JSR leak) or;
    * Run: Coupl_rep_uncoupl_Mult_clusters_mult_runs_save_1Dyad_StochJSR.m (randomly occurring JSR leak)
    * With same dependencies as **1** above
2. Result analysis & figure generation
    * Run: extract_Spont_spark_data.m (identifies and extracts spark statistics into one file)
    * Run: plot_Spont_Spark_Stats.m (swarm plots of spark amp and FDHM, statistical analysis)
    * With same dependencies as **2** above
