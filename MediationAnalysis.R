#---------------------------------------------------------------------------------------------------------------------
# All libraries here
library(boot)
library(survival)
library(data.table)
library(foreign)
library(dummies)
library(GenABEL)
library(dummies)
#---------------------------------------------------------------------------------------------------------------------
# Sources import here

# this script should be run from the same folder where src.R is
source('src.R')
#---------------------------------------------------------------------------------------------------------------------
# Define your parameters here!!!

#Data pathway
data_path<-"//storage.erasmusmc.nl/m/MyDocs/592004/My Documents/Desktop/Rscript4way/Test.sav"

#Path to save results
output<-'//storage.erasmusmc.nl/m/MyDocs/592004/My Documents/Desktop/Rscript4way/Test_results.csv'

  
#Define variables
A<<-'A2'
M<<-'M1'
Y<<-'Y1'
COVAR<<-c('C1','C2','C3')
  

  
#1=binary 0=continuous  
outcome=0
mediator=0

#Assign levels for the exposure that are being compared; 
#for mstar it is the level at which to compute the CDE and the remainder of the decomposition  
a<<-1     
astar<<-0 
mstar<<-0 

#Boostrap number of iterations
N_r=5

#---------------------------------------------------------------------------------------------------------------------
####### DONT TOUCH FROM HERE
#######
#---------------------------------------------------------------------------------------------------------------------

# Reading data file
data<-read.spss(data_path, to.data.frame=T) #TODO spss/csv/txt (?)


if (! prod(c(A,Y,M,COVAR) %in% names(data) ) )  {stop('Some of defined variable names are not in data file!')}

if ( mediator==1 & outcome==1 ) {  save_results(output=output, boot_function=boot.bMbO, N=N_r)  }
if ( mediator==0 & outcome==1 ) {  save_results(output=output, boot_function=boot.cMbO, N=N_r)  }
if ( mediator==1 & outcome==0 ) {  save_results(output=output, boot_function=boot.bMcO, N=N_r)  }
if ( mediator==0 & outcome==0 ) {  save_results(output=output, boot_function=boot.cMcO, N=N_r)  }




