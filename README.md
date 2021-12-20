# BayDERS
BAYesian Detection and Estimation of Rating Shifts

Purpose:
Software BayDERS (BAYesian Detection and Estimation of Rating Shifts) provides tools for the retrospective detection and estimation of stage-discharge rating shifts. These tools mainly use hydrometric information such as gaugings and stage record as input data to provide rating shift times with quantitative uncertainty. The user interface is available under R-studio. Source codes are written in R.


Acknowledegments:
The uncertainty quantification at mostly all steps of the tools in BayDERS is performed through BaM (BAyesian Modelling, developed by Benjamin Renard, INRAE with the use of a few fortran libraries developed by Dmitri Kavetski), which represents a probabilistic Bayesian approach and is based on the Monte Carlo Markov Chain algorithm (MCMC) to explore the posterior distribution of parameters. For a better description of this algorithm, please see Renard et al. [2006]. https://github.com/BaM-tools/BaM
The RC estimation is performed by means of BaRatin method, developed by Le Coz et al., 2014 (INRAE). https://riverhydraulics.inrae.fr/outils/logiciels/baratin-2/
The multi-period RC estimation is performed by means of BaRatin-SPD method, developed by Mansanarez et al., 2019 (INRAE). https://forge.irstea.fr/projects/bam/files


Disclaimer:
Please, notice that BayDERS is an experimental software. Further analysis is required for validating the proposed tools. We are not responsible for any loss (of data, profits, business, customers, orders, other costs or disturbances) derived by their use in the operational practice. The authorized user accepts the risks associated with using the software given its nature of free software. It is reserved to expert Users (developers or professionals) having prior knowledge in computer science, hydraulics and statistics. 

For any question or feedback please contact us at:
matteo.darienzo@cimafoundation.org, jerome.lecoz@inrae.fr or benjamin.renard@inrae.fr.
