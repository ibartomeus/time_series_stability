# time_series_stability

This repo explores the idea of assessing if a short term time series is compatible with a long term stable trend. We do that by using null models.

In folder `/R` we will need

1.  `get_params.R`: Function `get_params()` uses a inputs two vectors defining a time series. time and abundance. It extracts from this time serie its duration, variability, temporal autocorrelation and slope (save also model fit? Probably yes). Output as list of params / vector.

2.  `run_sims.R`: Function `run_sims()` takes as input time series parameters (duration, variability, temporal autocorrelation and slope) and creates `n` simulated time series with this parameters. Output vector of slopes and SE.

3.  `compare.R`: Function `compare()` that takes as in input the simulations using `run_sims()`, the observed value (using `get_params.R)`and calculates the Z-scores and p-values. Question, should we propagate error somehow?

In folder `/data` we will need:

1.  `BT_data.csv`: BIOIME2.0 time series with \> 20 years (restrict to insects?).

In `/scripts` folder:

1.  `test.Rmd`: This script should read one observed time series from BIOTIME. Calculate its parameters with `get_params()`, and run 100 simulations with slope = 0 and years = `n`. Select `n` random years (e.g. n = 6) and `get_params()` of this shorter time wondow. `compare()` the 100 slope simulations to the obtained value using Z-scores. Second, create a loop to go through all time series and store for each time series the original slope, the 6-years time window slope, its Z and p-values. Finally we are looking for a contingency table that for each combination of short term GLM and null model (rows: sig/ns, sig/sig, ns/sig, ns, ns) shows the number of time series with sig/ns Lont term trends (columns). If the method works, we expect for ST sig/ns to be LT ns, and for ST sig/sig to be LT sig. (ST ns/ns should be LTns, and STns/sig should be rare).
2.  Caveats: Using a random 6 years window might end up with low variability on responses. Using many time windows or selecting specific time windows (e.g. declining) might solve this.
3.  Further work: This might test that the null model approach works, BUT does not tell us anything about how to obtain the LT params when this info does not exists. This can be done using reference time-series, guesstimates, means for the target taxa elsewhere... To be discussed if the method works.
4.  Selling point: The key here is that we should interpret both tests toguether, the GLM and the Null model because they tell you different things. The first one talks about the ST trend, and the second one about how compatible is the observed ST trend with a LT stability.
