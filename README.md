# A 4-Way Decomposition
R code to implement the 4-way decomposition method on mediation and interaction for both binary and continuous mediator or outcome. This R code is easy-to-use by just giving a few inputs.

## User Guide
### Overview
The user must define the exposure, the mediator, the covariates, the outcome, and whether the mediator and outcome are binary or continuous. Also, the user must input the value of the exposure and mediator at which to compute the controlled direct effect and the remainder of the decomposition. Finally, the bootstrap number of iterations should be defined to obtain 95% confidence intervals. 
### Installation
`git clone` this repository (or save all files in separate folder).
### HOW TO USE
1.	Set variables in MediationAnalysis.R script **lines 15-43** according your analysis plan and save changes. 
2.	Run: 
3.	in the command line **Rscript** MediationAnalysis.R
4.	In Rstudio just run all code from MediationAnalysis.R
5.	copy/paste code from MediationAnalysis.R to R shell

We have set the mean value of the covariates at which the effects will be calculated. We considered the exposure/outcome and mediator/outcome to be confounded by the same set of confounders. We suggest the user to make a directed acyclic graph and check the assumptions for mediation analysis. 
## Citation
If you use this R code, please refer to the paper (including supplement) of VanderWeele (2014): ["A unification of mediation and interaction: a four-way decomposition."](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4220271/)
## License
This project is licensed under GNU GPL v3.

## Authors
* Unal Mutlu (Departments of Epidemiology and Ophthalmology, Erasmus MC, Rotterdam, Netherlands) 
* Gennady V. Roshchupkin (Department of Epidemiology, Radiology and Medical Informatics, Erasmus MC, Rotterdam, Netherlands)
* Gautam Sajeev (Department of Epidemiology, Harvard T.H. Chan School of Public Health, Boston, Massachusetts) 
* M. Kamran Ikram (Departments of Epidemiology and Neurology, Erasmus MC, Rotterdam, Netherlands)
* M. Arfan Ikram (Departments of Epidemiology, Neurology and Radiology, Erasmus MC, Rotterdam, Netherlands

## Contacts
If you have any questions/suggestions/comments or problems do not hesitate to contact us: u.mutlu@erasmusmc.nl, g.roshchupkin@erasmusmc.nl

