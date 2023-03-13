# Longitudinal-SST-analysis

Contains code for conducting Bayesian parametric estimation of Stop-signal task data in R. Requires the Dynamic Models of Choice package (https://osf.io/pbwx8/).

This analysis forms part of a larger project involving data from the Neuroimaging of the Children's Attention Project (NICAP; https://bmcpsychiatry.biomedcentral.com/articles/10.1186/s12888-016-0770-4).

The workflow is as follows: 

1. EXG3_DataPrep (Calculating traditional SSRT using ANALYZE-IT and prepping raw data for EXG3 Modelling)
2. EXG3_QA (QA and filterting of raw SST data)
3. EXG3_Modelling (Fitting QA'ed SST data to EXG3)
4. EXG3_MCMC (MCMC sampling of EXG3-fitted data to obtain subject-level posterior distributions for each parameter)
5. EXG3_Plotting (Plotting of posterior distributions)
6. EXG3_PostPred (Goodness-of-fit tests between the modelled and the observed data)

This analysis workflow employs the three-racer Ex-Gaussian Modelling paradigm deveopled by Prof. Dora Matzke and colleagues (2019) See link to the original publication: https://content.apa.org/record/2018-66881-003. Please ensure that you correctly cite all primary sources when attempting to use this code for your analysis.


For any questions pertaining to the use and dissemination of all scripts contained in this repository, please contact Mervyn Singh at mervynsingh87@gmail.com
