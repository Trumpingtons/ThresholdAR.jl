# Made by Pollis (Compose > New Methods > Package and Publish) from method.toml: change the method, not this file.
using Test, CSV, Downloads, Tables, ThresholdAR

const DATA = joinpath(@__DIR__, "data")

"An input table: a file of the tests, or downloaded from its link into memory."
function input(source; delim = ',', columns = nothing, header = true, skipto = nothing, limit = nothing)
	bytes = startswith(source, r"https?://") ? take!(Downloads.download(source, IOBuffer())) : read(joinpath(DATA, source))
	named = columns === nothing && header
	header = columns !== nothing ? Symbol.(columns) : named ? 1 : false
	rows = (; skipto = something(skipto, named ? 2 : 1), limit)
	file = delim == "whitespace" ? CSV.File(bytes; header, rows..., delim = ' ', ignorerepeated = true, stripwhitespace = true) : CSV.File(bytes; header, rows..., delim = delim == "tab" ? '\t' : only(delim))
	return Tables.columntable(file)
end

"The published table: its row labels and its values as printed, `nothing` where a cell is empty."
function target(file)
	table = CSV.File(joinpath(DATA, "targets", file * ".csv"); types = String, delim = ',')
	columns = Tables.columnnames(table)
	rows = [string(coalesce(row[columns[1]], "")) for row in table]
	values = [let cell = row[column]; cell === missing ? nothing : strip(cell) end for row in table, column in columns[2:end]]
	return rows, string.(columns[2:end]), values
end

"Tests each printed value against the result rounded to the digits printed."
function check(result, file, tolerance)
	rows, columns, values = target(file)
	result = result isa AbstractVecOrMat ? result : Tables.matrix(result)
	result = reshape(result, size(result, 1), :)
	@test size(result) == size(values)
	for i in axes(values, 1), j in axes(values, 2)
		printed = values[i, j]
		(printed === nothing || isempty(printed)) && continue
		mantissa, exponent = occursin(r"[eE]", printed) ? split(printed, r"[eE]") : (printed, "0")
		digits = (occursin('.', mantissa) ? length(split(mantissa, '.')[2]) : 0) - parse(Int, exponent)
		published = parse(Float64, printed)
		computed = round(result[i, j]; digits)
		@testset "$(rows[i]), $(columns[j])" begin
			@test isapprox(computed, published; atol = (tolerance + 1e-6) * 10.0^-digits)
		end
	end
end

@testset "ThresholdAR" begin
	@testset "Eq. 9.1, p. 259" begin
		lynx = input("https://raw.githubusercontent.com/Trumpingtons/PollisDatasets.jl/main/data/lynxmackenzie.csv")
		check(ThresholdAR.eq91(lynx), "eq91", 0)
	end
	@testset "Eq. 9.2, p. 261" begin
		sunspots = input("https://raw.githubusercontent.com/wch/r-source/7bcd838b66763d5029b5bdec668e3436c6e38017/src/library/datasets/data/sunspot.year.R"; columns = ["sunspots", "end"], skipto = 2, limit = 256)
		check(ThresholdAR.eq92(sunspots), "eq92", 0)
	end
	@testset "Eq. 9.3, p. 262" begin
		sunspots = input("https://raw.githubusercontent.com/wch/r-source/7bcd838b66763d5029b5bdec668e3436c6e38017/src/library/datasets/data/sunspot.year.R"; columns = ["sunspots", "end"], skipto = 2, limit = 256)
		check(ThresholdAR.eq93(sunspots), "eq93", 0)
	end
	@testset "Root mean RSS of Fig. 8 (p. 260) and RMSEs of Figs. 10 and 12 (pp. 261-262)" begin
		lynx = input("https://raw.githubusercontent.com/Trumpingtons/PollisDatasets.jl/main/data/lynxmackenzie.csv")
		sunspots = input("https://raw.githubusercontent.com/wch/r-source/7bcd838b66763d5029b5bdec668e3436c6e38017/src/library/datasets/data/sunspot.year.R"; columns = ["sunspots", "end"], skipto = 2, limit = 256)
		check(ThresholdAR.errors(lynx, sunspots), "errors", 1)
	end
	@testset "Lynx limit cycle, p. 259" begin
		lynx = input("https://raw.githubusercontent.com/Trumpingtons/PollisDatasets.jl/main/data/lynxmackenzie.csv")
		check(ThresholdAR.lynxcycle(lynx), "lynxcycle", 1)
	end
	@testset "Sunspot eventual forecasting function, p. 261" begin
		sunspots = input("https://raw.githubusercontent.com/wch/r-source/7bcd838b66763d5029b5bdec668e3436c6e38017/src/library/datasets/data/sunspot.year.R"; columns = ["sunspots", "end"], skipto = 2, limit = 256)
		check(ThresholdAR.sunspotcycle(sunspots), "sunspotcycle", 0)
	end
	@testset "Variances under eqs. 9.1-9.3 (pp. 259, 261, 262) and the NAICs of Figs. 8, 10 and 12 (pp. 260-262)" begin
		lynx = input("https://raw.githubusercontent.com/Trumpingtons/PollisDatasets.jl/main/data/lynxmackenzie.csv")
		sunspots = input("https://raw.githubusercontent.com/wch/r-source/7bcd838b66763d5029b5bdec668e3436c6e38017/src/library/datasets/data/sunspot.year.R"; columns = ["sunspots", "end"], skipto = 2, limit = 256)
		check(ThresholdAR.variances(lynx, sunspots), "variances", 1)
	end
end
