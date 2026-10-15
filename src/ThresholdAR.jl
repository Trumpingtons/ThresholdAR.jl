module ThresholdAR

using LinearAlgebra, Tables

export SETAR, setar, identify, selectdelay, residualvariances, predict, eventual, limitcycle, risedescent

"""
	SETAR

A self-exciting threshold autoregressive model with two regimes, SETAR(2; k₁, k₂), fitted by least squares:

	xₜ = a₀⁽ʲ⁾ + a₁⁽ʲ⁾ xₜ₋₁ + … + aₖⱼ⁽ʲ⁾ xₜ₋ₖⱼ + εₜ⁽ʲ⁾,   j = 1 if xₜ₋d ≤ r, j = 2 otherwise.

Fields: `d` (the delay), `r` (the threshold), `orders` (k₁, k₂), `coefficients` (a vector for each regime, the constant first),
`rss` and `n` (the residual sum of squares and the number of observations of each regime), `fitted` (the indices of the
observations fitted), `aic` (Tong and Lim's eq. 8.4, the sum over the regimes of nⱼ ln(RSSⱼ/nⱼ) + 2(kⱼ + 1)) and `naic` (the
AIC divided by the number of observations fitted, eq. 8.6).
"""
struct SETAR
	d::Int
	r::Float64
	orders::NTuple{2, Int}
	coefficients::NTuple{2, Vector{Float64}}
	rss::NTuple{2, Float64}
	n::NTuple{2, Int}
	fitted::UnitRange{Int}
	aic::Float64
	naic::Float64
end

Base.show(io::IO, m::SETAR) = print(io, "SETAR(2; $(m.orders[1]), $(m.orders[2])), d = $(m.d), r = $(m.r), NAIC = $(round(m.naic; digits = 5))")

"The regressors of an autoregression of order `k` at the times `times`: a constant, then xₜ₋₁ … xₜ₋ₖ."
design(x, times, k) = hcat(ones(length(times)), [x[t - i] for t in times, i in 1:k])

"The least-squares fit of an autoregression of order `k` at the times `times`, with its AIC (eq. 8.3)."
function regime(x, times, k)
	n = length(times)
	n > k + 1 || throw(ArgumentError("a regime has $n observations, too few for order $k"))
	A = design(x, times, k)
	b = A \ x[times]
	rss = sum(abs2, x[times] - A * b)
	(; k, b, rss, n, aic = n * log(rss / n) + 2(k + 1))
end

"The order from 0 to `maxorder` with the smallest AIC (eq. 8.3) and its fit."
function bestregime(x, times, maxorder)
	fits = [regime(x, times, k) for k in 0:min(maxorder, length(times) - 2)]
	fits[argmin([f.aic for f in fits])]
end

"The times fitted, from max(d, `maxorder`) + 1 to the end, and those of each regime."
function regimes(x, d, r, maxorder)
	fitted = max(d, maxorder) + 1:length(x)
	fitted, [t for t in fitted if x[t - d] <= r], [t for t in fitted if x[t - d] > r]
end

function model(d, r, fitted, f1, f2)
	aic = f1.aic + f2.aic
	SETAR(d, r, (f1.k, f2.k), (f1.b, f2.b), (f1.rss, f2.rss), (f1.n, f2.n), fitted, aic, aic / length(fitted))
end

"""
	setar(x, d, r, orders; maxorder = maximum(orders))

Fits SETAR(2; k₁, k₂) to the series `x` with delay `d`, threshold `r` and `orders` = (k₁, k₂) by least squares in each regime.
The fit starts at observation max(d, `maxorder`) + 1, so that models of orders up to `maxorder` are fitted to the same
observations and their AICs can be compared; give the largest order considered in the identification.
"""
function setar(x::AbstractVector{<:Real}, d::Integer, r::Real, orders; maxorder::Integer = maximum(orders))
	x = float.(x)
	fitted, times1, times2 = regimes(x, d, r, maxorder)
	model(d, r, fitted, regime(x, times1, orders[1]), regime(x, times2, orders[2]))
end

"""
	identify(x, d; maxorder, quantiles = 0.3:0.1:0.7)

Tong and Lim's identification of a SETAR(2; k₁, k₂) with delay `d` (Section 8): for each candidate threshold, the
`quantiles` of the threshold variable xₜ₋d over the observations fitted (the sorted value at position ⌊qN⌋), the order of
each regime is the one from 0 to `maxorder` with the smallest AIC (eq. 8.3); the threshold chosen is the one whose two
regimes have the smallest total AIC (eqs. 8.4 and 8.5).
"""
function identify(x::AbstractVector{<:Real}, d::Integer; maxorder::Integer, quantiles = 0.3:0.1:0.7)
	x = float.(x)
	fitted = max(d, maxorder) + 1:length(x)
	z = sort(x[fitted .- d])
	thresholds = unique([z[floor(Int, q * length(z))] for q in quantiles])
	models = map(thresholds) do r
		fitted, times1, times2 = regimes(x, d, r, maxorder)
		model(d, r, fitted, bestregime(x, times1, maxorder), bestregime(x, times2, maxorder))
	end
	models[argmin([m.aic for m in models])]
end

"""
	selectdelay(x; maxorder, delays = 1:maxorder, quantiles = 0.3:0.1:0.7)

The model identified (see [`identify`](@ref)) for each delay in `delays`, and the one with the smallest normalised AIC
(eq. 8.6), which compares fits to different numbers of observations. Returns `(; best, models)`. Tong and Lim also weigh
the eventual forecasting function and the one-step-ahead prediction errors of the observations held out before choosing.
"""
function selectdelay(x::AbstractVector{<:Real}; maxorder::Integer, delays = 1:maxorder, quantiles = 0.3:0.1:0.7)
	models = [identify(x, d; maxorder, quantiles) for d in delays]
	(; best = models[argmin([m.naic for m in models])], models)
end

"""
	residualvariances(m::SETAR)

The residual variance of each regime, RSSⱼ/nⱼ, and the pooled mean sum of squares of residuals, (RSS₁ + RSS₂)/N.
"""
residualvariances(m::SETAR) = [m.rss[1] / m.n[1], m.rss[2] / m.n[2], sum(m.rss) / sum(m.n)]

"The model's prediction of xₜ from x₁ … xₜ₋₁."
function step(m::SETAR, x, t)
	b = m.coefficients[x[t - m.d] <= m.r ? 1 : 2]
	b[1] + sum(b[i + 1] * x[t - i] for i in 1:length(b) - 1; init = 0.0)
end

"""
	predict(m::SETAR, x, times)

The one-step-ahead predictions of `x` at the indices `times`, each from the observations before it.
"""
predict(m::SETAR, x::AbstractVector{<:Real}, times) = [step(m, x, t) for t in times]

"""
	eventual(m::SETAR, x, steps)

The eventual forecasting function: the series `x` continued for `steps` steps by the model without noise.
"""
function eventual(m::SETAR, x::AbstractVector{<:Real}, steps::Integer)
	y = float.(collect(x))
	for _ in 1:steps
		push!(y, step(m, y, length(y) + 1))
	end
	y[length(x) + 1:end]
end

"""
	limitcycle(m::SETAR, x; steps = 5000, maxperiod = 150, tol = 1e-6)

The limit cycle the eventual forecasting function from `x` settles into after `steps` steps (Tong and Lim's Definition
7.1): one period of it, from its smallest value, or `nothing` if it has no period up to `maxperiod` (a limit point is a
cycle of period 1).
"""
function limitcycle(m::SETAR, x::AbstractVector{<:Real}; steps::Integer = 5000, maxperiod::Integer = 150, tol::Real = 1e-6)
	y = eventual(m, x, steps)
	period = findfirst(p -> maximum(abs.(y[end - p + 1:end] - y[end - 2p + 1:end - p])) < tol, 1:maxperiod)
	period === nothing && return nothing
	cycle = y[end - period + 1:end]
	circshift(cycle, 1 - argmin(cycle))
end

"""
	risedescent(cycle)

The local cycles of a limit cycle: for each local minimum, the rise to the next local maximum and the descent to the
following local minimum, in steps, as a vector of (rise, descent) pairs.
"""
function risedescent(cycle::AbstractVector{<:Real})
	p = length(cycle)
	c(i) = cycle[mod1(i, p)]
	minima = [i for i in 1:p if c(i - 1) > c(i) < c(i + 1)]
	maxima = [i for i in 1:p if c(i - 1) < c(i) > c(i + 1)]
	map(minima) do i
		top = minimum(j for j in [maxima; maxima .+ p] if j > i)
		bottom = minimum(j for j in [minima; minima .+ p; minima .+ 2p] if j > top)
		(top - i, bottom - top)
	end
end

# The tests: Tong and Lim's (1980) Section 9, on the Canadian lynx trappings (lynxmackenzie of PollisDatasets) and Wolf's
# yearly sunspot numbers (R's sunspot.year, the values the paper used)

lynxlog(lynx) = log10.(Float64.(lynx.lynx))
sunspots1700(sunspots, last) = Float64.(sunspots.sunspots[1:last - 1699])

lynx91(lynx) = setar(lynxlog(lynx), 2, 3.1163, (8, 3); maxorder = 10)
sunspots92(sunspots) = setar(sunspots1700(sunspots, 1920), 3, 36.6, (4, 12); maxorder = 20)
sunspots93(sunspots) = setar(sunspots1700(sunspots, 1924)[138:end], 5, 47.4, (4, 2); maxorder = 10)

"The coefficients of both regimes, the constant first, one row each."
coefficients(m) = [m.coefficients[1]; m.coefficients[2]]

"Eq. 9.1 (p. 259): the lynx SETAR(2; 8, 3), d = 2, fitted to 1831-1934."
eq91(lynx) = coefficients(lynx91(lynx))

"Eq. 9.2 (p. 261): the sunspot SETAR(2; 4, 12), d = 3, fitted to 1720-1920."
eq92(sunspots) = coefficients(sunspots92(sunspots))

"Eq. 9.3 (p. 262): the sunspot SETAR(2; 4, 2), d = 5, fitted to 1847-1924."
eq93(sunspots) = coefficients(sunspots93(sunspots))

"""
The lynx limit cycle (p. 259), in the order printed: the paper lists it backwards from its smallest value (read forwards,
it rises for 3 years and falls for 6, but the text gives a rise of six years and a descent of three).
"""
function lynxcycle(lynx)
	cycle = limitcycle(lynx91(lynx), lynxlog(lynx))
	[cycle[1]; reverse(cycle[2:end])]
end

"The sunspot eventual forecasting function (p. 261): its period, then the rise and descent of each local cycle."
function sunspotcycle(sunspots)
	cycle = limitcycle(sunspots92(sunspots), sunspots1700(sunspots, 1920))
	[length(cycle); collect(Iterators.flatten(risedescent(cycle)))]
end

"""
The variances printed under eqs. 9.1-9.3 and the NAICs of Figs. 8, 10 and 12. They are those of the same models with one
lag fewer in each regime, and the NAIC penalises them as if they had the orders stated: the coefficients are not.
"""
function shortfit(x, d, r, orders, maxorder)
	m = setar(x, d, r, orders .- 1; maxorder)
	[residualvariances(m); (m.aic + 2 * 2) / length(m.fitted)]
end

"Variances and NAIC printed with eq. 9.1 (p. 259) and Fig. 8 (p. 260)."
lynxprinted(lynx) = shortfit(lynxlog(lynx), 2, 3.1163, (8, 3), 10)

"Variances and NAIC printed with eq. 9.2 (p. 261) and Fig. 10 (p. 261)."
sunspots92printed(sunspots) = shortfit(sunspots1700(sunspots, 1920), 3, 36.6, (4, 12), 20)

"Variances and NAIC printed with eq. 9.3 (p. 262) and Fig. 12 (p. 262)."
sunspots93printed(sunspots) = shortfit(sunspots1700(sunspots, 1924)[138:end], 5, 47.4, (4, 2), 10)

"The variances and NAICs printed with eqs. 9.1, 9.2 and 9.3, one column each."
variances(lynx, sunspots) = [lynxprinted(lynx) sunspots92printed(sunspots) sunspots93printed(sunspots)]

rmse(e) = sqrt(sum(abs2, e) / length(e))

"""
The fit statistics of Figs. 8, 10 and 12: the lynx root mean RSS (1831-1934), and the root-mean-square errors of the
one-step-ahead sunspot predictions of 1921-1955 (eq. 9.2) and 1925-1944 (eq. 9.3).
"""
function errors(lynx, sunspots)
	m = lynx91(lynx)
	x = Float64.(sunspots.sunspots)
	years(y) = y .- 1699
	[sqrt(sum(m.rss) / length(m.fitted)),
		rmse(x[years(1921:1955)] - predict(sunspots92(sunspots), x, years(1921:1955))),
		rmse(x[years(1925:1944)] - predict(sunspots93(sunspots), x, years(1925:1944)))]
end

end
