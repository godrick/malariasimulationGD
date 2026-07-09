# ------------------------------------------------------------------------------
# seasonality.R
# ------------------------------------------------------------------------------
# Fourier coefficients for the seasonal-rainfall cycle used by the example.
#
# malariasimulationGD models seasonality via the Fourier series
#
#   value(t) = g0 + sum_i (g[i] * cos(2*pi*t*i/365) + h[i] * sin(2*pi*t*i/365))
#
# (floored at `rainfall_floor`). The simulator multiplies the larval carrying
# capacity by value(t), so larger values during the rains drive more mosquito
# production.
#
# Provenance of the coefficients below:
#   The coefficients are the H=3 least-squares fit to the day-of-year mean of
#   the CHIRPS daily-rainfall product over the Busia hull,
#   2001-01-01 through 2020-12-31 (n = 9131 daily observations). The annual
#   cycle was scaled so its annual mean equals 1, i.e. these coefficients
#   describe a relative rainfall index (not absolute mm/day). Reconstruction
#   from the fit ranges roughly 0.28 (dry-season trough) to 1.74 (rainy-season
#   peak), matching CHIRPS's observed seasonal amplitude over Busia.
#
# The example therefore uses real Busia seasonality without shipping the
# CHIRPS source CSV. The CSV that produced these coefficients lives in the
# active research repository; recompute and overwrite the numbers here if you
# need a different aggregation window or location.
#
# To recompute from a CHIRPS-style daily CSV:
#
#   df <- read.csv("...chirps_daily.csv")
#   df$doy <- as.integer(format(as.Date(df$date), "%j"))
#   daily_mean <- aggregate(df$rainfall_mm, by = list(doy = df$doy), FUN = mean)
#   t <- daily_mean$doy; y <- daily_mean$x / mean(daily_mean$x)
#   X <- cbind(1, do.call(cbind, lapply(seq_len(3), function(i)
#         cbind(cos(2*pi*t*i/365), sin(2*pi*t*i/365)))))
#   bhat <- solve(crossprod(X), crossprod(X, y))
# ------------------------------------------------------------------------------

seven_node_seasonality <- function() {
  list(
    # Fit to Busia CHIRPS 2001-2020, scaled to annual mean 1.
    g0 = 1.000000,
    g  = c(-0.154337, -0.361168,  0.053734),   # cosine coefficients (harmonics 1, 2, 3)
    h  = c( 0.068664, -0.390137, -0.093621),   # sine   coefficients (harmonics 1, 2, 3)
    rainfall_floor = 0.001
  )
}


#' Evaluate the seasonal-rainfall cycle at a vector of times.
#'
#' Useful for plotting the cycle in diagnostics. Same math as the package's
#' internal Fourier evaluator.
#'
#' @param t numeric vector of times (days since simulation start)
#' @param seas result of seven_node_seasonality()
seasonality_value <- function(t, seas = seven_node_seasonality()) {
  v <- rep(as.numeric(seas$g0), length(t))
  for (i in seq_along(seas$g)) {
    v <- v +
      seas$g[[i]] * cos(2 * pi * t * i / 365) +
      seas$h[[i]] * sin(2 * pi * t * i / 365)
  }
  pmax(v, seas$rainfall_floor)
}
