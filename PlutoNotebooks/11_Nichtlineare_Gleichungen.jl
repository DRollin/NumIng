### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ afca956e-5c1f-11ed-1721-815ee5a6146c
begin
	using PlutoUI, Plots, ColorSchemes, ForwardDiff, LaTeXStrings

	H(x) = x < 0 ? 0 : 1
	precompile(H, (Function, Float64))
	precompile(H, (Function, Float32))
	precompile(H, (Function, Float16))

	const arccos = acos
	const arcsin = asin
	const arctan = atan
	
	md"""
	# Nicht-Lineare Gleichungen

	Hinweis:

	Für die Darstellungen werden die Betrachteten Funktionen für Argumente außerhalb des Definitionsbereiches der Funktion zu ``0`` gestezt.
	"""
end

# ╔═╡ ee3773b6-7d0e-4d4a-93fc-7f50bce0fec7
begin
	maxiter = 20
	x = -10:0.1:10
	
	md"""
	## Lösen einer Gleichung in Nullstellen-Form

	Betrachtete Gleichung: ``f(x) =`` $(@bind fstring TextField(; default="cos(x) + x^2")) ``= 0`` für ``x \in [-10,10]``

	Zum Lösen der Gleichung werden das Newton- und das Sekanten-Verfahren verwendet mit:
	
	- Anzahl Iterationen: ``N_{iter} =`` $(@bind Nᵢₜₑᵣ Slider(1:1:maxiter; default=1, show_value=true))
	- Startwert: ``x^{(0)}`` = $(@bind x₀ Slider(x; default=9, show_value=true))
	- Sekanten-Verfahren: ``x^{(-1)}`` = $(@bind x₋₁ Slider(x; default=10, show_value=true))
	"""
end

# ╔═╡ 0ea08f90-f347-496d-8491-2302a6a424f0
begin
	f = eval(Meta.parse("x -> x ≥ x[1] && x ≤ x[end] ? " * fstring * " : 0.0"))
	∂ₓf(x) = ForwardDiff.derivative(f, x)
	
	xₙ = zeros(maxiter+1)
	xₙ[1] = x₀
	
	xₛ = zeros(maxiter+2)
	xₛ[1] = x₋₁
	xₛ[2] = x₀
	
	for n in 1:maxiter
		xₙ[n+1] = xₙ[n] - f(xₙ[n]) / ∂ₓf(xₙ[n])
		n += 1
		xₛ[n+1] = xₛ[n] - f(xₛ[n]) *(xₛ[n] - xₛ[n-1]) / (f(xₛ[n]) - f(xₛ[n-1]))
	end

	md"""
	Die Iterationsvorschriften der Verfahren lauten:

	- Newton: ``x^{(k+1)} = x^{(k)} - \frac{ f(x^{(k)}) }{ f'(x^{(k)}) }``
	- Sekanten: ``x^{(k+1)} = x^{(k)} - f(x^{(k)}) \frac{ x^{(k)} - x^{(k-1)} }{ f(x^{(k)}) - f(x^{(k-1)}) }``
	"""
end

# ╔═╡ 416aafca-472b-4e18-b4ab-4386bbadf1b0
begin
	pₙ = plot(x, f.(x); label="f(x)", title="Newton-Verfahren", framestyle=:zerolines, legend=:topleft, linewidth=2, xlabel=L"x", ylabel=L"f")
	scatter!(pₙ, xₙ[1:Nᵢₜₑᵣ+1], f.(xₙ[1:Nᵢₜₑᵣ+1]); label=nothing)
	for n in 1:Nᵢₜₑᵣ
		plot!(pₙ, [xₙ[n], xₙ[n+1]], [f(xₙ[n]), 0]; color=:green, arrow=true, label=false)
		plot!(pₙ, [xₙ[n+1], xₙ[n+1]], [0, f(xₙ[n+1])]; color=:grey, arrow=true, label=false)
	end

	pₛ = plot(x, f.(x); label=L"f(x)", title="Sekanten-Verfahren", framestyle=:zerolines, legend=:topleft, linewidth=2, xlabel=L"x", ylabel=L"f")
	scatter!(pₛ, xₛ[1:Nᵢₜₑᵣ+2], f.(xₛ[1:Nᵢₜₑᵣ+2]); label=false)
	for n in 2:Nᵢₜₑᵣ+1
		plot!(pₛ, [xₛ[n-1], xₛ[n], xₛ[n+1]], [f(xₛ[n-1]), f(xₛ[n]), 0]; color=:green, arrow=true, label=false)
		plot!(pₛ, [xₛ[n+1], xₛ[n+1]], [0, f(xₛ[n+1])]; color=:grey, arrow=true, label=false)
	end
	md"""
	$(plot(pₙ, pₛ; layout=(2,1), size=(700, 400*2)))
	"""
end

# ╔═╡ b69ae320-214f-4ae7-9976-5abf3e7384d8
begin
	maxiterᵖ = 20
	xᵖ = -10:0.1:10
	md"""
	## Lösen einer Gleichung in Fixpunkt-Form

	Betrachtete Fixpunkt-Form: ``T(x) =`` $(@bind Tstring TextField(;default="atan(x) +(x/3)^2-4")) ``= x`` für ``x \in [-10,10]``

	Zum Lösen der Gleichung wird das Fixpunkt-Verfahren verwendet mit:
	
	- Anzahl Iterationen: ``N_{iter} =`` $(@bind Nᵢₜₑᵣᵖ Slider(1:1:maxiterᵖ; default=1, show_value=true))
	- Startwert: ``x^{(0)} =`` $(@bind x₀ᵖ Slider(-10:0.1:10; default=-9, show_value=true))
	"""
end

# ╔═╡ 11e00be8-9264-4560-a127-e1919536d6e5
begin
	T = eval(Meta.parse("x -> x ≥ xᵖ[1] && x ≤ xᵖ[end]  ? " * Tstring * " : 0.0"))

	xₚ = zeros(maxiterᵖ+1)
	xₚ[1] = x₀ᵖ
	for n in 1:maxiterᵖ
		xₚ[n+1] = T(xₚ[n])
	end

	md"""
	Iteriert wird nach der Vorschrift **xₙ₊₁ = T(xₙ)**
	"""
end

# ╔═╡ 92df17cc-1d7b-4bb6-91cc-044d523c8cb6
begin
	pₚ = plot(xᵖ, T.(xᵖ); label=L"T(x)", title="Fixpunkt-Iteration", framestyle=:zerolines, legend=:topleft, linewidth=2, size=(700,400), xlabel=L"x", ylabel=L"T")
	plot!(pₚ, xᵖ, xᵖ; color=:grey, label=L"x")
	scatter!(pₚ, xₚ[1:Nᵢₜₑᵣᵖ+1], T.(xₚ[1:Nᵢₜₑᵣᵖ+1]); label=false, color=:orange)
	
	for n in 1:Nᵢₜₑᵣᵖ
		plot!(pₚ, [xₚ[n], xₚ[n+1]], [T(xₚ[n]), xₚ[n+1]]; color=:green, arrow=true, label=false)
		plot!(pₚ, [xₚ[n+1], xₚ[n+1]], [xₚ[n+1], T(xₚ[n+1])]; color=:grey, arrow=true, label=false)
	end
	
	pₚ
end

# ╔═╡ Cell order:
# ╟─afca956e-5c1f-11ed-1721-815ee5a6146c
# ╟─ee3773b6-7d0e-4d4a-93fc-7f50bce0fec7
# ╟─0ea08f90-f347-496d-8491-2302a6a424f0
# ╟─416aafca-472b-4e18-b4ab-4386bbadf1b0
# ╟─b69ae320-214f-4ae7-9976-5abf3e7384d8
# ╟─11e00be8-9264-4560-a127-e1919536d6e5
# ╟─92df17cc-1d7b-4bb6-91cc-044d523c8cb6
