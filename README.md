# ThresholdAR.jl

Threshold autoregression: SETAR models and limit cycles.

Tong and Lim's self-exciting threshold autoregressive model with two regimes, SETAR(2; k1, k2): least-squares fit, identification of the threshold and the orders by AIC, choice of the delay by normalised AIC, one-step-ahead prediction, and the eventual forecasting function with its limit cycle.

Article: Tong, H. and Lim, K. S. (1980). Threshold autoregression, limit cycles and cyclical data. Journal of the Royal Statistical Society: Series B 42(3), 245-268. [doi:10.1111/j.2517-6161.1980.tb01126.x](https://doi.org/10.1111/j.2517-6161.1980.tb01126.x)

## Usage

```julia
using ThresholdAR, CSV, Tables
lynx = Tables.columntable(CSV.File(download("https://raw.githubusercontent.com/Trumpingtons/PollisDatasets.jl/main/data/lynxmackenzie.csv")))
x = log10.(Float64.(lynx.lynx))  # the Canadian lynx trappings of 1821-1934, logarithms to the base 10
selectdelay(x; maxorder = 10, delays = 1:6).models  # the model identified for each delay, with its NAIC
m = setar(x, 2, 3.1163, (8, 3); maxorder = 10)  # Tong and Lim's eq. 9.1
m.coefficients, residualvariances(m)
cycle = limitcycle(m, x)  # the limit cycle of the eventual forecasting function
risedescent(cycle)  # the rise and descent of each local cycle
```

## Vignette

[Threshold autoregression: SETAR models and limit cycles](vignette/ThresholdAR.pdf)

## Tests

The tests compute the published tables below from their data and check every printed value, rounded to the digits printed:

- Eq. 9.1, p. 259 (`test/data/targets/eq91.csv`)
- Eq. 9.2, p. 261 (`test/data/targets/eq92.csv`)
- Eq. 9.3, p. 262 (`test/data/targets/eq93.csv`)
- Root mean RSS of Fig. 8 (p. 260) and RMSEs of Figs. 10 and 12 (pp. 261-262) (`test/data/targets/errors.csv`)
- Lynx limit cycle, p. 259 (`test/data/targets/lynxcycle.csv`)
- Sunspot eventual forecasting function, p. 261 (`test/data/targets/sunspotcycle.csv`)
- Variances under eqs. 9.1-9.3 (pp. 259, 261, 262) and the NAICs of Figs. 8, 10 and 12 (pp. 260-262) (`test/data/targets/variances.csv`)

The data, downloaded from its links when the tests run or kept in `test/data/inputs`:

- `lynx`: <https://raw.githubusercontent.com/Trumpingtons/PollisDatasets.jl/main/data/lynxmackenzie.csv>, PollisDatasets: the number of lynx trapped in the Mackenzie River district of Canada each year from 1821 to 1934 (Elton and Nicholson, 1942)
- `sunspots`: <https://raw.githubusercontent.com/wch/r-source/7bcd838b66763d5029b5bdec668e3436c6e38017/src/library/datasets/data/sunspot.year.R>, R's sunspot.year: Wolf's yearly sunspot numbers from 1700, the values Tong and Lim used (SILSO has since revised them); the first 256, 1700 to 1955. Fetched, not kept: SILSO's licence (CC BY-NC 4.0) does not allow redistribution with Pollis

Run them with `using Pkg; Pkg.test("ThresholdAR")`.

Made with [Pollis](https://pollis.org).
