#Rscript CIV_full_analysis_V1.R ${model} ${sim} ${gen} $exp_idx ${mtCN} ${mutRate}
#model: linear_bn - bottleneckmodel, linear_const - constant model
#sim: simulation ID
#gen: generation [50, 100, 150, 200, 250, 300, 350, 400]
#exp_idx: index for the intensity of clonal expansions. 1: 0.1, 
#                                                       2: 0.2,
#                                                       3: 0.25,
#                                                       4: 0.5,
#                                                       5: 0.75,
#                                                       6: 0.8,
#                                                       7: 0.9
#mtCN: copy number of mtDNA [500, 750, 1000] 
#mutRate: mutation rate of mtDNA [1e-08, 5e-08, 1e-07]


#example run
Rscript CIV_full_analysis_V1.R linear_bn 257859 400 1 500 5e-08
