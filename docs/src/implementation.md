# Implementation and results

In this page, we solve numerically the optimal control problem presented in the [home page](index.md) for an environment with n = 3 substrates.

```@example diauxic
using DifferentialEquations
using Plots
using OptimalControl
using NLPModelsIpopt
using LaTeXStrings
using Plots.PlotMeasures
```

## Problem definition

The rate functions are defined as Michaelis–Menten kinetics:

```math
w_i(s_i) = k \frac{s_i}{K_i + s_i}, \quad i = 1, 2, 3, \qquad w_R(m) = k \frac{m}{K_R + m},
```

where ``K_i`` is the half-saturation constant of substrate uptake (common to the three substrates), ``K_R`` is the half-saturation constant of protein synthesis, and ``k = 1`` in both examples below, so that the functions saturate at 1. Non-negativity constraints on the states are also added to help the solver.

```@example diauxic
# Define uptake functions
wm1(s) = k*s/(Kᵢ + s)
wm2(s) = k*s/(Kᵢ + s)
wm3(s) = k*s/(Kᵢ + s)
wR(m) = k * m / (KR + m)

function diauxic_ocp(t₀, tf, I₀, Y₁, Y₂, Y₃; scale = 1, tol = 1e-12)
    ocp = @def begin
        # time interval
        t ∈ [t₀, tf], time
        # state variables: z = (s₁, s₂, s₃, m, x)  
        z = (s₁, s₂, s₃, m, x) ∈ R⁵, state   
        # control variables: u = (u₀, u₁, u₂, u₃)
        u = (u₀, u₁, u₂, u₃) ∈ R⁴, control   
        # Control bounds: each u_i ∈ [0,1]
        [0, 0, 0, 0] ≤ u(t) ≤ [1, 1, 1, 1]    
        u₀(t) + u₁(t) + u₂(t) + u₃(t) ≤ 1 # Control bounds: sum of controls ≤ 1

        z(t) ≥ [0, 0, 0, 0, 0]           # State bounds: all states non-negative

        z(t₀) == [I₀[1], I₀[2], I₀[3], I₀[4], I₀[5]]  # Initial conditions

        ż(t) == [
            - u₁(t) * wm1(s₁(t)) * x(t),                          
            - u₂(t) * wm2(s₂(t)) * x(t),                            
            - u₃(t) * wm3(s₃(t)) * x(t),                          
              Y₁*u₁(t)*wm1(s₁(t)) + Y₂*u₂(t)*wm2(s₂(t))         
              + Y₃*u₃(t)*wm3(s₃(t)) - u₀(t)*wR(m(t))*(m(t)+1),                                 
              u₀(t)*wR(m(t)) * x(t)            
        ]
        scale * x(tf) → max     # Maximize final biomass
    end 
    
    sol = solve(
        ocp, :collocation, :adnlp, :ipopt;
        scheme = :gauss_legendre_3,
        grid_size   = 1000,       
        tol         = tol,
        display     = false
        )

    return sol 
end
nothing # hide
```

## First example: substrate-abundant case

In the substrate-abundant case, the substrates in the medium are initially abundant with respect to the intermediate metabolites in the cells, so that investing in metabolism is initially more profitable than investing in growth. We consider three substrates with decreasing yields Y₁ > Y₂ > Y₃ and initial conditions satisfying ``Y_1 w_1(s_{10}) > Y_2 w_2(s_{20}) > Y_3 w_3(s_{30})``, so that the substrate with the highest yield is the most profitable one at the beginning of the process. The initial amount of intermediate metabolites is small, so the cell is expected to first invest all its resources in metabolism. This corresponds to the first general trajectory presented in [Yabo (2025)](https://doi.org/10.1016/j.automatica.2024.111881).

```@example diauxic
# Initial conditions and time horizon
s₁₀, s₂₀, s₃₀, m₀, x₀ = 1e-3, 2e-3, 3e-3, 1e-3, 5e-3 
t₀, tf = 0.0, 2.0
Y₁, Y₂, Y₃ = 1.0, 0.5, 0.2
# Define uptake parameters
k, KR, Kᵢ = 1, 0.003, 0.0003
I₀ = [s₁₀, s₂₀, s₃₀, m₀, x₀]

# Solve the optimal control problem
sol = diauxic_ocp(t₀, tf, I₀, Y₁, Y₂, Y₃)
nothing # hide
```

```@example diauxic
sol # hide
```

### Switching times

According to the theoretical analysis of [Yabo (2025)](https://doi.org/10.1016/j.automatica.2024.111881), the optimal solution is a concatenation of three types of arcs, characterized by the value of u₀:

- **pure-metabolism arcs** (u₀ = 0), where all resources are allocated to the metabolic machinery and biomass is not produced;
- **pure-growth arcs** (u₀ = 1), where all resources are allocated to ribosomes, and substrate uptake is arrested;
- **singular arcs** (0 < u₀ < 1), where resources are shared between metabolism and growth.

Moreover, along the singular arc, the *i*-th substrate starts being consumed when ``Y_i w_i(s_i)`` reaches the maximum of these functions over the substrates already being consumed. These switching times are detected numerically from the optimal solution and displayed in the plots below.

```@setup diauxic
function calculate_switching_times(sol, Y₁, Y₂, Y₃)
  
    arc1, arc2, arcr1, arcr2 = nothing, nothing, nothing, nothing
    
    try 
        t = time_grid(sol)
        u = control(sol)
    
        # Evaluate controls at time points
        u0_vals = [u(t)[1] for t in t]
        
        # Beginning and end of the singular arc
        if Y₁ ≈ 1.0 && Y₂ ≈ 0.3 && Y₃ ≈ 0.1
            arcr1_idx = findfirst(u0_vals .< 0.9)
            arcr2_idx = findlast(u0_vals .< 0.9)
        else
            arcr1_idx = findfirst(u0_vals .> 0.1)
            arcr2_idx = findfirst(u0_vals .> 0.5)
        end
        
        # Activation of u₂ and u₃: first time each control exceeds a small threshold
        ε = 0.05
        arc1_idx = findfirst(ti -> u(ti)[3] > ε, t)
        arc2_idx = findfirst(ti -> u(ti)[4] > ε, t)

        if arc1_idx !== nothing
            arc1 = t[arc1_idx]
        end
        if arc2_idx !== nothing
            arc2 = t[arc2_idx]
        end
        
        if arcr1_idx !== nothing
            arcr1 = t[arcr1_idx]
            println("arcr1 =", arcr1)
        end
        if arcr2_idx !== nothing
            arcr2 = t[arcr2_idx]
            println("arcr2 =", arcr2)
        end
    catch e
        println("Error in switching time calculation: ", e)
    end
    
    return arc1, arc2, arcr1, arcr2
end

# Calculate switching times for our solution
arc1, arc2, arcr1, arcr2 = calculate_switching_times(sol, Y₁, Y₂, Y₃)
```

### States

We now plot the optimal state trajectories. In all plots, the shaded regions at the beginning and at the end of the time interval correspond to the non-singular arcs, and the dashed vertical lines indicate the instants at which s₂ and s₃ start being consumed.

```@example diauxic
function plot_state_variables(sol, t₀, tf, arc1, arc2, arcr1, arcr2)    
    # Plot 1: Substrates
    substrates_plot = plot(
        t -> state(sol)(t)[3], 0, tf,
        label = L"$s_3$", color = "#79af97", lw = 5,
        grid = true,
        gridlinewidth = 1.5,
        size = (1300, 700),
        legendfontsize = 22,
        tickfontsize = 20,
        guidefontsize = 23,
        framestyle = :box,
        left_margin = 10mm,
        bottom_margin = 10mm,
        foreground_color_subplot = :black, 
        xlabel = "t",
        ylabel = "Extracellular quantities"
    )
    plot!(substrates_plot,
        t -> state(sol)(t)[2], 0, tf,
        label = L"$s_2$", color = "#df8f44", lw = 5
    )
    plot!(substrates_plot,
        t -> state(sol)(t)[1], 0, tf,
        label = L"$s_1$", color = "#b24746", lw = 5
    )
    vspan!(substrates_plot, [t₀, arcr1], fillalpha=0.2, fillcolor=:grey, label="")
    vspan!(substrates_plot, [arcr2, tf], fillalpha=0.2, fillcolor=:grey, label="")
    vline!(substrates_plot, [arc1], color=:grey, linestyle=:dash, linewidth=5, label="")
    vline!(substrates_plot, [arc2], color=:grey, linestyle=:dash, linewidth=5, label="")

    # Plot 2: Metabolites
    metabolites_plot = plot(
        t -> state(sol)(t)[4], 0, tf,
        label = L"$m$", color = "#374e55", lw = 5,
               grid = true,
        gridlinewidth = 1.5,
        size = (1300, 700),
        legendfontsize = 22,
        tickfontsize = 20,
        guidefontsize = 23,
        framestyle = :box,
        left_margin = 10mm,
        bottom_margin = 10mm,
        foreground_color_subplot = :black, 
        xlabel = "t",
        ylabel = "metabolite"
    )
    vline!(metabolites_plot, [arc1], color=:grey, linestyle=:dash, linewidth=5, label="")
    vline!(metabolites_plot, [arc2], color=:grey, linestyle=:dash, linewidth=5, label="")
    vspan!(metabolites_plot, [t₀, arcr1], fillalpha=0.2, fillcolor=:grey, label="")
    vspan!(metabolites_plot, [arcr2, tf], fillalpha=0.2, fillcolor=:grey, label="")
    # Plot 3: Biomass
    biomass_plot = plot(
        t -> state(sol)(t)[5], 0, tf,
        label = L"$x$", color = "#374e55", lw = 5,
        grid = true,
        gridlinewidth = 1.5,
        size = (1300, 700),
        legendfontsize = 22,
        tickfontsize = 20,
        guidefontsize = 23,
        framestyle = :box,
        left_margin = 10mm,
        bottom_margin = 10mm,
        foreground_color_subplot = :black, 
        xlabel = "t",
        yformatter = :plain,
        ylabel = "biomass"
    )

    vline!(biomass_plot, [arc1], color=:grey, linestyle=:dash, linewidth=5, label="")
    vline!(biomass_plot, [arc2], color=:grey, linestyle=:dash, linewidth=5, label="")
    vspan!(biomass_plot, [t₀, arcr1], fillalpha=0.2, fillcolor=:grey, label="")
    vspan!(biomass_plot, [arcr2, tf], fillalpha=0.2, fillcolor=:grey, label="")
    return substrates_plot, metabolites_plot, biomass_plot
end

# Generate the state variable plots
substrates_plot, metabolites_plot, biomass_plot = plot_state_variables(sol, t₀, tf, arc1, arc2, arcr1, arcr2)
nothing # hide
```

```@example diauxic
substrates_plot
```

The substrates are consumed sequentially, starting from the one with the highest yield, which is the diauxic behaviour described in the [home page](index.md).

```@example diauxic
metabolites_plot
```

The intermediate metabolites first accumulate during the initial pure-metabolism arc, are kept at a regulated level along the singular arc, and are finally converted into biomass during the final pure-growth arc.

```@example diauxic
biomass_plot
```

### Controls

The optimal controls are plotted stacked on top of each other, so that the total height of the colored area is u₀ + u₁ + u₂ + u₃ = 1, as all the available resources are used along the optimal solution. The fraction allocated to ribosomes u₀ is zero along the initial pure-metabolism arc, takes intermediate values along the singular arc, and is equal to one along the final pure-growth arc. Meanwhile, the metabolic controls u₁, u₂ and u₃ are activated one after the other, at the switching times computed above.

```@example diauxic
function plot_controls(sol, t₀, tf; arc1=nothing, arc2=nothing, arcr1=arcr1, arcr2=arcr2)
    # Extract time grid and control function
    t_vals = time_grid(sol)
    u = control(sol)
    
    # # Evaluate controls at time points
    u0_vals = [u(t)[1] for t in t_vals]
    u1_vals = [u(t)[2] for t in t_vals]
    u2_vals = [u(t)[3] for t in t_vals] 
    u3_vals = [u(t)[4] for t in t_vals]
     
    p = plot(t_vals, u0_vals, 
         fillrange=0, fillalpha=0.6, fillcolor="#374e55",
         label=L"$u_0$", linewidth=0, 
         grid=true, gridlinewidth=1.5,
         size=(1600, 1000),
         legendfontsize=24, tickfontsize=20, guidefontsize=26,
         framestyle=:box, left_margin=10mm, bottom_margin=10mm,
         foreground_color_subplot=:black,
         xlabel="t", ylabel="Optimal controls",
         xlims=(t₀, tf), ylims=(0, 1)
    )
    
    # u1 on top of u0
    plot!(p, t_vals, u0_vals .+ u1_vals, 
          fillrange=u0_vals, fillalpha=0.9, fillcolor="#b24746", 
          label=L"$u_1$", linewidth=0)
    
    # u2 on top of u0+u1  
    plot!(p, t_vals, u0_vals .+ u1_vals .+ u2_vals, 
          fillrange=(u0_vals .+ u1_vals), fillalpha=0.9, fillcolor="#df8f44", 
          label=L"$u_2$", linewidth=0)
    
    # u3 on top of u0+u1+u2
    plot!(p, t_vals, u0_vals .+ u1_vals .+ u2_vals .+ u3_vals, 
          fillrange=(u0_vals .+ u1_vals .+ u2_vals), fillalpha=0.9, 
          fillcolor="#79af97", label=L"$u_3$", linewidth=0)
    
    plot!(p, t_vals, u0_vals, linewidth=0.5, linecolor=:black, alpha=1, label="")
    plot!(p, t_vals, u0_vals .+ u1_vals, linewidth=0.5, linecolor=:black, label="")
    plot!(p, t_vals, u0_vals .+ u1_vals .+ u2_vals, linewidth=0.5, linecolor=:black, label="")
    
    # Adding vertical lines for switching times
    if arc1 !== nothing
        vline!(p, [arc1], color=:grey, linestyle=:dash, linewidth=5, alpha=0.9, label="")
    end
    if arc2 !== nothing
        vline!(p, [arc2], color=:grey, linestyle=:dash, linewidth=5, alpha=0.9, label="")
    end
    if arcr1 !== nothing
        vline!(p, [arcr1], color=:grey, linestyle=:dash, linewidth=5, alpha=0.9, label="")
    end
    if arcr2 !== nothing
        vline!(p, [arcr2], color=:grey, linestyle=:dash, linewidth=5, alpha=0.9, label="")
    end

    plot!(p, legend=:topright, legend_columns=2)
    
    return p
end

# Plot the controls for our solution
control_plot = plot_controls(sol, t₀, tf; arc1=arc1, arc2=arc2, arcr1=arcr1, arcr2=arcr2)
```

## Second example: metabolite-abundant case

In the metabolite-abundant case, the situation is reversed: the intermediate metabolites are initially abundant with respect to the substrates, so that investing in growth is initially more profitable than investing in metabolism. We now consider three substrates with the same initial concentration and a high initial amount of intermediate metabolites m₀. In this situation, the cell is expected to first invest all its resources in growth, in order to convert the available metabolites into biomass, which produces an initial pure-growth arc. A larger pool of metabolites also increases the duration of the final pure-growth arc. This corresponds to the second general trajectory presented in [Yabo (2025)](https://doi.org/10.1016/j.automatica.2024.111881). Since the substrate concentrations are an order of magnitude smaller than in the first example, the objective is very flat with respect to the activation times of the controls; we thus scale the objective to improve the numerical conditioning of the problem.

```@example diauxic
k, KR, Kᵢ = 1, 0.3, 0.015
s₁₀ = s₂₀ = s₃₀ = 1e-4
Y₁, Y₂, Y₃ = 1.0, 0.3, 0.1
m₀, x₀ = 0.1, 0.02
t₀, tf = 0, 5

I₀_new = [s₁₀, s₂₀, s₃₀, m₀, x₀]

# Solve with new parameters
sol2 = diauxic_ocp(t₀, tf, I₀_new, Y₁, Y₂, Y₃; scale = 1e4)
nothing # hide
```

```@example diauxic
sol2 # hide
```

```@setup diauxic
# Calculate switching times for the new parameter set
arc1_new, arc2_new, arcr1_new, arcr2_new = calculate_switching_times(sol2, Y₁, Y₂, Y₃)
```

```@example diauxic
substrates_plot_bis, metabolites_plot_bis, biomass_plot_bis = plot_state_variables(sol2, t₀, tf, arc1_new, arc2_new, arcr1_new, arcr2_new)
nothing # hide
```

```@example diauxic
substrates_plot_bis
```

```@example diauxic
metabolites_plot_bis
```

```@example diauxic
biomass_plot_bis
```
