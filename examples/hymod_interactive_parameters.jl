### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 71c490a0-506a-11eb-0f26-b3924a24895f
using HTTP, CSV, DataFrames, Hymod, Plots

# ╔═╡ ce0c5324-506b-11eb-1855-159926367ef0
begin
	gr()
	theme(:bright)
end

# ╔═╡ 83363da2-506a-11eb-1bfd-23e7c9ea4ab9

df = begin
	# define url to test data
	testDataUrl = "https://raw.githubusercontent.com/JuliaHydro/Hymod.jl/master/test/data/test_forcings.csv"
	# get response
	# read response as a DataFrame
	CSV.File(HTTP.get(testDataUrl).body) |> DataFrame
end

# ╔═╡ 92b9384c-506a-11eb-0904-554f3ca6fe86
# create a column in the data frame for potential evapotranspiration
df[:,:pet] = hargreaves(df,tminCol=:tmin,tmaxCol=:tmax,dtCol=:date)

# ╔═╡ a4d4c2a0-506d-11eb-3937-75e927332299
@bind bexp html"<input type='range' min='0' max='2' value='0.9' step='0.1'> <label>bexp</label>"

# ╔═╡ d86d17de-506d-11eb-0fd1-c7acae2ced4a
@bind kq html"<input type='range' min='0.5' max='1.2' value='0.9' step='0.1'> <label>kq</label>"

# ╔═╡ d990b6fc-506d-11eb-0815-2bd333ce5ae8
@bind cmax html"<input type='range' min='0.01' max='1' value='0.9' step='0.01'> <label>cmax</label>"

# ╔═╡ 204f1464-506c-11eb-3f6e-6df9b9b0db43
@bind alpha html"<input type='range' min='0.2' max='0.99' value='0.9' step='0.01'> <label>alpha</label>"

# ╔═╡ daae38ac-506d-11eb-317d-4dd9e5319d51
@bind ks html"<input type='range' min='0.01' max='0.5' value='0.4' step='0.01'> <label>ks</label>"

# ╔═╡ e6793368-506c-11eb-1d23-5d97cf757030
pars = Dict(
	:alpha => alpha,
	:bexp => bexp,
	:kq => kq,
	:cmax => cmax,
	:ks =>ks
)

# ╔═╡ a1b5f8bc-506a-11eb-12b3-e9f89e677320
# run a simulation
q = simulate(df,precipCol=:precip, petCol=:pet; pars...)

# ╔═╡ aca2edac-506a-11eb-01e0-6bb9898a1879
begin
	plot(q, label="Simulated", dpi=450)
	plot!(df[:,:obs], label="Observed")
end

# ╔═╡ Cell order:
# ╠═71c490a0-506a-11eb-0f26-b3924a24895f
# ╠═ce0c5324-506b-11eb-1855-159926367ef0
# ╠═83363da2-506a-11eb-1bfd-23e7c9ea4ab9
# ╠═92b9384c-506a-11eb-0904-554f3ca6fe86
# ╠═e6793368-506c-11eb-1d23-5d97cf757030
# ╠═a1b5f8bc-506a-11eb-12b3-e9f89e677320
# ╠═a4d4c2a0-506d-11eb-3937-75e927332299
# ╠═d86d17de-506d-11eb-0fd1-c7acae2ced4a
# ╠═d990b6fc-506d-11eb-0815-2bd333ce5ae8
# ╠═204f1464-506c-11eb-3f6e-6df9b9b0db43
# ╠═daae38ac-506d-11eb-317d-4dd9e5319d51
# ╠═aca2edac-506a-11eb-01e0-6bb9898a1879
