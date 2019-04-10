##################################################

Paper: 
A Flexible Replication-Based Classification Approach for Parkinson's Disease Detection by Using Voice Recordings

Authors:
Lizbeth Naranjo (1), Ruth Fuentes-Garcia (1), Carlos J. Perez (2).

(1) Universidad Nacional Autonoma de Mexico, Mexico
(2) Universidad de Extremadura, Spain

Title: 
33rd National Forum of Statistics (FNE) and 13th Latin-American Congress of Statistical Societies (CLATSE)

Publisher:
Springer Proceedings in Mathematics & Statistics


##################################################

Instructions to run the codes in R and JAGS are provided. 
The codes are applied to obtain a similar analysis as in Section 5 ‘Results’, but without cross-validation. 

##################################################

##################################################
FILES 

The file 'ParkinsonMixtureR.R' contains the R code. The JAGS code is run from this R file.

The file 'ParkinsonMixtureJAGS.bug' contains the JAGS model. 

The used dataset can be downloaded from
UCI Machine Learning Repository (https://archive.ics.uci.edu/ml/datasets/NUESTRA_BASE_PARKINSON**).

##################################################

To run the files, do the following.
 
1.- Download JAGS from www.mcmc-jags.sourceforge.net/

2.- Install the packages necessary to run the R file. 
These are indicated in the R file. 

3.- Change the address indicated in ‘setwd()’. 
setwd("HERE"). This is the address where the file ‘ParkinsonMixtureJAGS.bug’ is in.

4.- Read the data. The data file should be located in the indicated address.


##################################################
Related references (using the same dataset)

- Naranjo,L., Perez,C.J., Campos-Roca,Y., Martın,J. (2016).
 Addressing voice recording replications for Parkinson's disease detection. 

Expert Systems With Applications 46, 286–292

- Naranjo, L.,Perez,C.J., Martın,J., Campos-Roca,Y. (2017).
 A two-stage variable selection and classification approach for Parkinson's disease detection by using voice recording replications. 
Computer Methods and Programs in Biomedicine 142, 147-156.

##################################################


