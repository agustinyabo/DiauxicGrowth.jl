# Application to experimental data

In this page, we calibrate the model presented in the [home page](index.md) with experimental data of *E. coli* growing on glucose and lactose[^1], following the case study of [^2]. The model has n = 2 substrates, where s₁ and s₂ represent glucose and lactose, respectively. Instead of solving an optimal control problem, the model is simulated in closed loop with a feedback control law derived from the optimal solutions, and its parameters are fitted to the measurements by least squares.

```@example calibration
using DifferentialEquations
using Optim
using Plots
using LaTeXStrings
using Plots.PlotMeasures
```

## Experimental data

The data consists of 15 measurements of the concentrations of glucose, lactose and biomass, in g/L, along a batch process of about 6.5 hours.

```@example calibration
# Time [h], glucose, lactose and biomass concentrations [g/L]
t_data   = [0.0, 0.5, 1.0, 1.5, 1.9166, 2.6166, 3.05, 3.5, 3.9166, 4.25, 4.7166, 5.0833, 5.45, 5.95, 6.5833]
glc_data = [0.21741594111, 0.19652882013, 0.16994521161, 0.13196862801, 0.08924497146, 0.00759531672, 0.0,
            0.00094941459, 0.00094941459, 0.00189882918, 0.0, -0.00284824377, 0.00189882918, 0.00094941459, 0.00379765836]
lac_data = [1.140143985, 1.1681803125, 1.170049401, 1.1644421355, 1.1532276045, 1.1307985425, 1.080333153,
            0.9924859935, 0.8840788605, 0.721468161, 0.52334478, 0.314006868, 0.041119947, 0.0, -0.0056072655]
x_data   = [0.0315909, 0.0402066, 0.0561616, 0.0749885, 0.0938154, 0.1483815, 0.1678466, 0.1956083,
            0.2370913, 0.2772979, 0.3369696, 0.3832391, 0.4295086, 0.4735444, 0.4898185]
nothing # hide
```

The measurements show the typical diauxic pattern: glucose is consumed first, while lactose remains constant, and lactose only starts being consumed once glucose is exhausted, with a change in the growth rate of the biomass at the switch.

```@example calibration
style = (grid = true, gridlinewidth = 1.5, legendfontsize = 16, tickfontsize = 14, guidefontsize = 17,
         framestyle = :box, foreground_color_subplot = :black)

data_plot = plot(t_data, x_data; color = "#374e55", lw = 2, marker = :circle, ms = 7, markeralpha = 0.6, label = "Biomass",
                    xlabel = "t [h]", ylabel = "Concentrations [g/L]", legend = :left,
                    size = (1200, 700), left_margin = 10mm, bottom_margin = 10mm, style...)
plot!(data_plot, t_data, glc_data; color = "#b24746", lw = 2, marker = :circle, ms = 7, markeralpha = 0.6, label = "Glucose")
plot!(data_plot, t_data, lac_data; color = "#79af97", lw = 2, marker = :circle, ms = 7, markeralpha = 0.6, label = "Lactose")
```

The states of the model are non-dimensional: concentrations are converted to the model variables by multiplying them by the parameter β = 0.003 L/g, which represents the inverse of the cell density. The initial conditions of the substrates and the biomass are taken from the first measurement.

```@example calibration
const β, K_R = 0.003, 1.0
z₀(m₀) = [glc_data[1] * β, lac_data[1] * β, m₀, x_data[1] * β]
nothing # hide
```

## Feedback control law

As the final pure-growth arc of the optimal solutions becomes negligible for long processes, the optimal controls can be approximated by a feedback law, obtained by neglecting the final arc of the optimal control structure[^2].

!!! tip "Theorem"
    Let ``\rho(s) = \max(Y_1 w_1(s_1), \dots, Y_n w_n(s_n))`` and ``\mathcal{I}(s) = \left\{ i = 1, \dots, n \, : \, \rho(s) = Y_i w_i(s_i) \right\}``. The feedback control law is

    ```math
    \begin{aligned}
    u_0(s, m, x) &= \begin{cases}
    0 & \text{if } \rho(s) > \dfrac{w_R^2(m)}{w_R'(m)}, \\[2mm]
    1 & \text{if } \rho(s) < \dfrac{w_R^2(m)}{w_R'(m)}, \\[2mm]
    u_{0,\rm sing}(s, m, x) & \text{if } \rho(s) = \dfrac{w_R^2(m)}{w_R'(m)},
    \end{cases} \\[4mm]
    u_i(s, m, x) &= \begin{cases}
    \dfrac{1 - u_0(s, m, x)}{\sum_{j \in \mathcal{I}(s)} \frac{w_i'(s_i)}{w_j'(s_j)}} & \text{if } i \in \mathcal{I}(s), \\[2mm]
    0 & \text{otherwise},
    \end{cases}
    \end{aligned}
    ```

    where the singular control is

    ```math
    u_{0,\rm sing}(s, m, x) = \frac{x + \phi(s, m) \frac{w_R(m)}{w_R'(m)}}{x + \phi(s, m) \left(m + 1 + \frac{w_R(m)}{w_R'(m)}\right)},
    ```

    with

    ```math
    \phi(s, m) = \left(2 w_R'(m) - \frac{w_R(m)}{w_R'(m)} w_R''(m)\right) \sum_{j \in \mathcal{I}(s)} \frac{1}{w_j'(s_j)}.
    ```

In our case, the metabolic machinery is allocated to glucose only, until ``Y_1 w_1(s_1)`` reaches ``Y_2 w_2(s_2)``. From then on, both substrates are consumed simultaneously. In the implementation, the equalities defining the singular arc and the set ``\mathcal{I}(s)`` are checked up to a small tolerance. The rate functions are Michaelis–Menten kinetics with half-saturation constants ``K_1`` and ``K_2`` for glucose and lactose, and ``K_R`` for protein synthesis. Time is rescaled by the maximal rate ``k``, so that the time in hours is t/k.

```@example calibration
function feedback(z, θ; singular_only = false)
    _, K₁, K₂, Y₁, Y₂ = θ
    s₁, s₂, m, x = z
    K, Ks₁, Ks₂ = β * K_R, β * K₁, β * K₂
    w₁, w₂, wR = s₁ / (Ks₁ + s₁), s₂ / (Ks₂ + s₂), m / (K + m)
    w₁′, w₂′ = (1 - w₁) / (Ks₁ + s₁), (1 - w₂) / (Ks₂ + s₂)
    wR′ = (1 - wR) / (K + m)
    wR″ = -2wR′ / (K + m)

    # Distribution of the metabolic machinery between glucose and lactose
    δ₁ = (Y₁ * w₁ - Y₂ * w₂)^2 < 1e-4 ? w₂′ / (w₁′ + w₂′) : 1.0
    δ₂ = 1 - δ₁

    # Fraction of resources allocated to ribosomes
    maxYw, threshold = max(Y₁ * w₁, Y₂ * w₂), wR^2 / wR′
    if singular_only || (maxYw - threshold)^2 < 1e-3
        a = 2wR′ - wR / wR′ * wR″
        b = w₁′ * δ₁ * x
        u₀ = (b + a * wR / wR′) / (b + a * (m + 1 + wR / wR′))
    else
        u₀ = maxYw > threshold ? 0.0 : 1.0
    end
    return u₀, δ₁ * (1 - u₀), δ₂ * (1 - u₀)
end

function dynamics!(dz, z, p, t)
    θ, singular_only = p
    _, K₁, K₂, Y₁, Y₂ = θ
    s₁, s₂, m, x = z
    w₁, w₂, wR = s₁ / (β * K₁ + s₁), s₂ / (β * K₂ + s₂), m / (β * K_R + m)
    u₀, u₁, u₂ = feedback(z, θ; singular_only)
    dz[1] = -u₁ * w₁ * x
    dz[2] = -u₂ * w₂ * x
    dz[3] = Y₁ * u₁ * w₁ + Y₂ * u₂ * w₂ - u₀ * wR * (m + 1)
    dz[4] = u₀ * wR * x
end

function simulate(θ; m₀ = 0.001, singular_only = false, kwargs...)
    prob = ODEProblem(dynamics!, z₀(m₀), (0.0, θ[1] * t_data[end]), (θ, singular_only))
    return solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-12, maxiters = 100_000, verbose = false, kwargs...)
end
nothing # hide
```

## Calibration

The parameters θ = (k, K₁, K₂, Y₁, Y₂) are adjusted by minimizing the squared difference between the model and the measurements:

```math
f(\theta) = \sum_{i=1}^{15} \left[ \left(\frac{s_1(t_i)}{\beta} - \hat s_1^i\right)^2 + \left(\frac{s_2(t_i)}{\beta} - \hat s_2^i\right)^2 + \left(\frac{x(t_i)}{\beta} - \hat x^i\right)^2 \right],
```

where ``\hat s_1^i``, ``\hat s_2^i`` and ``\hat x^i`` are the *i*-th measurements of glucose, lactose and biomass. The bang arcs of the feedback law make this cost function discontinuous with respect to the parameters, which prevents the optimizer from converging. For the calibration, we thus use a smooth version of the law, in which u₀ is always given by the singular control, with an initial amount of intermediate metabolites m₀ = 0.01.

```@example calibration
function cost(θ; kwargs...)
    sol = simulate(θ; saveat = θ[1] .* t_data, kwargs...)
    SciMLBase.successful_retcode(sol) || return Inf
    Z = reduce(hcat, sol.u)
    return sum(abs2, Z[1, :] / β .- glc_data) +
           sum(abs2, Z[2, :] / β .- lac_data) +
           sum(abs2, Z[4, :] / β .- x_data)
end

calibration_cost(θ) = cost(θ; m₀ = 0.01, singular_only = true)

lower = [1.0, 0.01, 0.01, 0.1, 0.1]   # k, K₁, K₂, Y₁, Y₂
upper = [10.0, 2.0, 2.0, 1.0, 1.0]
θ_init = [2.0, 0.2, 0.2, 0.9, 0.2]

result = optimize(calibration_cost, lower, upper, θ_init, Fminbox(LBFGS()))
θ_fit = Optim.minimizer(result)
```

The calibrated parameters are close to the values reported in [^2] (k = 1.77 h⁻¹, K₁ = 0.037 g/L, K₂ = 0.01 g/L, Y₁ = 0.77 and Y₂ = 0.25).

## Comparison with the experimental data

We now simulate the model with the calibrated parameters and the complete feedback law, including the initial bang arc, with m₀ = 0.001.

```@example calibration
sol = simulate(θ_fit; saveat = 0.02)
k = θ_fit[1]
t = sol.t ./ k                         # time in hours
Z = reduce(hcat, sol.u)
U = reduce(hcat, [collect(feedback(z, θ_fit)) for z in sol.u])

# End of the initial pure-metabolism arc and activation of the lactose uptake
t_arc0 = t[findfirst(U[1, :] .> 0.001)]
t_lac = t[findfirst(U[3, :] .> 0.05)]
nothing # hide
```

The top plot compares the simulation with the measurements, while the bottom plot shows the controls of the feedback law. The shaded region indicates the initial pure-metabolism arc, and the dashed line the instant at which lactose starts being consumed.

```@example calibration
p1 = plot(t, Z[4, :] / β; color = "#374e55", lw = 4, label = L"x/\beta", ylabel = "Concentrations [g/L]",
          xlims = (0, t_data[end]), ylims = (0, 1.25), legend = :left, legend_columns = 2, style...)
scatter!(p1, t_data, x_data; color = "#374e55", alpha = 0.6, ms = 7, label = "Biomass")
plot!(p1, t, Z[1, :] / β; color = "#b24746", lw = 4, label = L"s_1/\beta")
scatter!(p1, t_data, glc_data; color = "#b24746", alpha = 0.6, ms = 7, label = "Glucose")
plot!(p1, t, Z[2, :] / β; color = "#79af97", lw = 4, label = L"s_2/\beta")
scatter!(p1, t_data, lac_data; color = "#79af97", alpha = 0.6, ms = 7, label = "Lactose")

p2 = plot(t, U[1, :]; fillrange = 0, fillalpha = 0.7, fillcolor = "#374e55", linewidth = 0, label = L"u_0",
          xlabel = "t [h]", ylabel = "Controls", xlims = (0, t_data[end]), ylims = (0, 1), legend = :right, style...)
plot!(p2, t, U[1, :] .+ U[2, :]; fillrange = U[1, :], fillalpha = 0.9, fillcolor = "#b24746", linewidth = 0, label = L"u_1")
plot!(p2, t, U[1, :] .+ U[2, :] .+ U[3, :]; fillrange = U[1, :] .+ U[2, :], fillalpha = 0.9, fillcolor = "#79af97",
      linewidth = 0, label = L"u_2")

for p in (p1, p2)
    vspan!(p, [0, t_arc0]; fillalpha = 0.3, fillcolor = :grey, label = "")
    vline!(p, [t_lac]; color = :grey, linestyle = :dash, linewidth = 3, label = "")
end

plot(p1, p2; layout = grid(2, 1, heights = [0.67, 0.33]), size = (1200, 1100), left_margin = 10mm, bottom_margin = 5mm)
```

The model reproduces the diauxic growth pattern: glucose, which has the highest yield, is consumed first, and lactose uptake starts when ``Y_1 w_1(s_1)`` reaches ``Y_2 w_2(s_2)``. The largest deviation from the data occurs around this switch: in the experiment, glucose is fully exhausted before lactose starts being consumed, while in the model both substrates are consumed simultaneously after the switch. As discussed in [^2], this is because the model does not allow for a complete depletion of the substrates in finite time.

[^1]: **Bettenbrock, K., Fischer, S., Kremling, A., Jahreis, K., Sauter, T., & Gilles, E. D.** (2006). A quantitative approach to catabolite repression in *Escherichia coli*. *Journal of Biological Chemistry*, 281(5), 2578-2584. [https://pubmed.ncbi.nlm.nih.gov/16263707/](https://pubmed.ncbi.nlm.nih.gov/16263707/)

[^2]: **Yabo, A. G.** (2025). Optimal control strategies in a generic class of bacterial growth models with multiple substrates. *Automatica*, 171, 111881. [https://doi.org/10.1016/j.automatica.2024.111881](https://doi.org/10.1016/j.automatica.2024.111881)
