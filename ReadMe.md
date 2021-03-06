## Management Strategy Evaluation using ADMB
Steve Martell 

International Pacific Halibut Commission

___

This is a small demonstration model for developing the infrastructure for conducting management strategy evaluation on a multicore unix based machine running AD Model Builder.

* * *

Has a cool Web-based interface!

![Alt text](./MSEE.png)

See the code documentation at http://smartell.github.io/MSEdemo/

#### Makefile targets

* `make`		   : compile's code and runs operating model as estimator
* `make mse`       : runs the operating model and runs MSE in the FINALS_SECTION
* `make clean`     : clean up temporary and model output files in the directory.
* `make data`      : create simulation directories and copy files to each directory.
* `make sims`      : runs MSE simulations (tip: use -j8 to run 8 jobs in parallel).
* `make collect`   : run R-script to produce allsims.Rdata for use in R.
* `make cleansims` : permanently remove simulation directories.