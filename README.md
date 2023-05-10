# tas
Code used to perform the analyses in Rohlfs, Chris, 2023. "Forbidden Knowledge and Specialized Training: A Versatile Solution for the Two Main Sources of Overfitting in Linear Regression," The American Statistician 77(2): 160-8. https://doi.org/10.1080/00031305.2022.2128874

baseline.empirical.R - Computes the baseline specification for the empirical application
	specifications:
		line 7: cores <- 20 #for multicore processing.
			Note that, because parallel processing was used, seed does not suffice to exactly replicate paper results.
		line 177: pct <- 0.50 #determines fraction of sample to be used for training
	inputs:
		mri.RDS - cleaned predictors and outcomes from the Neurocognitive Aging Data.
	outputs:
		components.RDS - large (~0.6 gig) file containing a data.table of results by observation, iteration, & spec
		rss.RDS - smaller data.table with summary statistics for each specification & iteration
		MRI_50pct_rss.png - sample output graph; more in doc.graphs.R

baseline.simulation.R - Computes the baseline specification for the Monte Carlo simulation
	specifications:
		lines 6 to 13 contain specification details.
		Note: for homoskedastic specification, data generation & specification details are slightly different. A sample version of that calculation appears in homosked45_125.R

	outputs:
		simulation.data.hetero.txt - text version of components (as in empirical). Saved as text so that "append" can be used given the large size of the file. Some versions (depending on parameters) get as large as 25-30 gigabytes

		sims.RDS - smaller data.table with summary statistics for each specification & iteration

homosked45_125.R - Computes an alternative homoskedastic specification for Monte Carlo
	- Note: also calculates some extra statistics related to MAPE (used in figures) and sigmas

doc.graphs - Generates all the graphs for the paper. Requires that the empirical & simulation scripts have been run with various parameter configurations.
