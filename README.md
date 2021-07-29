# pf - pointforecasts

This is a replication repository for the article "Interpretation of point forecats with unknown directive" by Patrick Schmidt, Matthias Katzfuss and Tilmann Gneiting.

For an accesible documentation of the method see the R-Package PointFore on [CRAN](https://cran.r-project.org/web/packages/PointFore/index.html) or [Github](https://github.com/Schmidtpk/PointFore).

This repository uses the PointFore Package to generate the simulations in the paper. The simulations are executed and results saved with the R-scripts:
* starndardR2_serverx.R (x=1,2)
* Rtwostep_server.R (x=1,2,3,4)
* specifications.R
* specificationsExpetiles.R

and the results are loaded and plotted in

* standardR_serverresults.R
* twostepR_serverresults.R
* specifications_serverresults.R

For questions please write to pschmidte@gmail.com.