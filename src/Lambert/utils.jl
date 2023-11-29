#============================================
minLambert: THIS IS BROKEN

Description: solve for the minimum total delta V solution by iterating through Lambert outputs

Inputs: 
    1. x0_P - initial state of the pursuer at t0, [DU; DU/TU]
    2. x0_E - initial state of evader at t0, [DU; DU/TU]
    3. mu - gravitational parameter in DU^3/TU^2
    4. timeLims - time limits of the optimization, [t0, tf] in TU

Outputs: 
    1. Minimum total delta V solution
============================================#
function minLambert(x0_P, x0_E, mu, timeLims)

    # Lambert solution finds the velocity vectors (v1, v2) that connects the 
    # position vectors (r1, r2)
    # This function iterates through different dm, and times (bounded by timeLims) to improve the 
    # (v1, v2) found by lambertbattin()

    # generate time vector to test
    tStart = LinRange(timeLims[1]+0.0001, timeLims[2], 100)
    tEnd = LinRange(timeLims[1]+0.0001, timeLims[2], 100)

    # test both prograde and retrograde solutions
    dm = ["pro", "retro"]
    
    for i in 1:length(dm)
        # stack output data for plotting
        output = []

        for j in 1:length(tStart)
            
            # propagate pursuer to start time
            xt_P_t, xt_P_u = propagate_2Body(x0_P, tStart[j], mu)
            xt_P = xt_P_u[end]
            rt_P = xt_P[1:3]
            vt_P = xt_P[4:6]

            # propagate evader to start time
            xt_E = propagate_2Body(x0_E, tStart[j], mu).u[end]

            for k in 1:length(tEnd)
                # only test if end time is greater than start time
                if tEnd[j] < tStart[j]
                    continue
                end

                # propagate evader to end time
                xf_E_t, xf_E_u = propagate_2Body(xt_E, tEnd[j], mu)
                xf_E = xf_E_U[end]

                # get final state components for evader
                rf_E = xf_E[1:3]
                vf_E = xf_E[4:6]

                # get lambert solution between initial pursuer state to final evader state
                v1, v2 = lambertbattin(rt_P, rf_E, mu, dm[i], tEnd[j])

                dv1 = norm(vt_P - v1)
                dv2 = norm(vf_E - v2)

                push!(output, [dv1, dv2, tStart[j], tEnd[j]])
            end
 
        end
        # plotDV(output, dm[i])
        
        # println(dm[i])
        # println("Δv1 = ", minimum(dv1_vals), " [km/s]")
        # println("Δv2 = ", minimum(dv2_vals), " [km/s]")
    end

    return output
end

export lambertContour
#============================================
plotDV:

Description: generate contour plot to vary both t0 and tf 

Inputs: 

Outputs: 
    1. Plot for prograde solutions, log(|Δv|)
    2. Plot for prograde solutions, |Δv|
    3. Plot for retrograde solutions, log(|Δv|)
    4. Plot for retrograde solutions, |Δv|
============================================#
function lambertContour(t, x₀_P, x₀_E, mu)

    t0 = 0.0001:t[2]/100:t[2]
    tf = 0.0001:t[2]/100:t[2]

    dm = ["pro", "retro"]
    output = []

    for i in 1:length(dm)
        dir = dm[i]
        vals = []
        for i=1:length(t0)
            for j=1:length(tf)
                if t0[i] >= tf[j]
                    continue
                end
        
                dv1, dv2 = testLambert(x₀_P, x₀_E, mu, t0[i], tf[j], dir)
                push!(vals, [t0[i], tf[j], dv1+dv2])
            end
        end
        
        vals = mapreduce( permutedims, vcat, vals) 
        
        fig = Figure()
        ax = Axis(fig[1, 1]; aspect=DataAspect())
        c = contourf!(ax, vals[:,1], vals[:,2], log.(vals[:,3]), levels = range(0, 100, length = 100), extendlow = :cyan, extendhigh = :magenta)
        Colorbar(fig[1, 2], c, label="log(|Δv|) (km/s)")
        ax.xlabel = "Maneuver 1 time (s)"
        ax.ylabel = "Maneuver 2 time (s)"
        if dir == "pro"
            ax.title = "Prograde Lambert Solution: log(|Δv|) (km/s)"
        else
            ax.title = "Retrograde Lambert Solution: log(|Δv|) (km/s)"
        end
        save("plots/lambert_"*dir*"_logdvTotal_contour.png", fig)
        
        fig = Figure()
        ax = Axis(fig[1, 1]; aspect=DataAspect())
        c = contourf!(ax, vals[:,1], vals[:,2], vals[:,3], levels = range(0, 100, length = 100), extendlow = :cyan, extendhigh = :magenta)
        Colorbar(fig[1, 2], c, label="|Δv| (km/s)")
        ax.xlabel = "Maneuver 1 time (s)"
        ax.ylabel = "Maneuver 2 time (s)"
        if dir == "pro"
            ax.title = "Prograde Lambert Solution: |Δv| (km/s)"
        else
            ax.title = "Retrograde Lambert Solution: |Δv| (km/s)"
        end
        
        save("plots/lambert_"*dir*"_dvTotalcontour.png", fig)

        push!(output, vals)
    end

    return output
end 

export varyT0
#============================================
varyT0: 

Description: Generate plot to show how dv1 and dv2 change as a function of t0, with a fixed tf

Inputs: 
    1. kep0_P - initial keplerian elements of pursuer at t0
    2. x0_E - initial state of evader at t0, [DU; DU/TU]
    3. t - time limits, [t0, tf] in TU
    4. mu - gravitational parameter in DU^3/TU^2

Outputs:
    1. Plot of trajectories generated by prograde solutions
    2. Plot of trajectories generated by retrograde solutions
    3. Plot of trajectories with both solutions
    4. Plot of tStart vs dv1
    5. Plot of tStart vs dv1, dv2
    6. Plot of tStart vs dv1 + dv2

============================================#
function varyT0(kep0_P, x₀_E, t, mu)

    # define times between t0 and tf
    tStart = LinRange(t[1], t[2], 100)
    tStart = tStart[2:end] # skip first time step to since propagator won't do anything

    # define initial pursuer state
    x0_P = kep2cart(kep0_P, mu)

    # propagate final evader state
    x_E_t, x_E_u = propagate_2Body(x₀_E, tStart[end], mu)
    x_E = x_E_u[end]
    rf_E = x_E[1:3]
    vf_E = x_E[4:6]

    # convert to matrix for plotting
    x_E = mapreduce(permutedims, vcat, x_E_u)

    dv1_pro = []
    dv1_retro = []
    dv2_pro = []
    dv2_retro = []

    traj_P = []
    traj_prograde = []
    traj_retrograde = []

    for j in 1:length(tStart)
        currentStart = tStart[j]

        # propagate pursuer to tStart
        x_P_t, x_P_u = propagate_2Body(x0_P, currentStart, mu)
        rt_P = x_P_u[end][1:3]
        vt_P = x_P_u[end][4:6]
        
        # get lambert solution between pursuer state at t to evader state at tStart[end]
        v1_pro, v2_pro = lambertbattin(rt_P, rf_E, mu, "pro", tStart[end])
        v1_retro, v2_retro = lambertbattin(rt_P, rf_E, mu, "retro", tStart[end])

        # compute dv magnitudes
        push!(dv1_pro, norm(vt_P - v1_pro))
        push!(dv1_retro, norm(vt_P - v1_retro))
        push!(dv2_pro, norm(vf_E - v2_pro))
        push!(dv2_retro, norm(vf_E - v2_retro))

        x_P = mapreduce(permutedims, vcat, x_P_u)
        push!(traj_P, x_P)

        # generate output for prograde trajectories        
        x₀_P_lambert = [rt_P; v1_pro]
        x_P_lambert_t, x_P_lambert_u = propagate_2Body(x₀_P_lambert, tStart[end], mu)
        x_P_lambert = mapreduce(permutedims, vcat, x_P_lambert_u) 
        push!(traj_prograde, x_P_lambert)

        # generate output for retrograde trajectories
        x₀_P_lambert = [rt_P; v1_retro]
        x_P_lambert_t, x_P_lambert_u = propagate_2Body(x₀_P_lambert, tStart[end], mu)
        x_P_lambert = mapreduce(permutedims, vcat, x_P_lambert_u) 
        push!(traj_retrograde, x_P_lambert)
    end

    generateLambertPlots(x_E, traj_P, traj_prograde, traj_retrograde, tStart, dv1_pro, dv1_retro, dv2_pro, dv2_retro)
    return [x_E, traj_P, traj_prograde, traj_retrograde, tStart, dv1_pro, dv1_retro, dv2_pro, dv2_retro]
end

export varyTF
#============================================
varyTF: 

Description: Generate plot to show how dv1 and dv2 change as a function of tf, with a fixed t0

Inputs: 
    1. kep0_P - initial keplerian elements of pursuer at t0
    2. x0_E - initial state of evader at t0, [DU; DU/TU]
    3. t - time limits, [t0, tf] in TU
    4. mu - gravitational parameter in DU^3/TU^2

Outputs:
    1. Plot of trajectories generated by prograde solutions
    2. Plot of trajectories generated by retrograde solutions
    3. Plot of trajectories with both solutions
    4. Plot of tStart vs dv1
    5. Plot of tStart vs dv1, dv2
    6. Plot of tStart vs dv1 + dv2

============================================#
function varyTF(kep0_P, x₀_E, t, mu)

    # define initial pursuer state
    x0_P = kep2cart(kep0_P, mu)

    # define times between t0 and tf
    tEnd = LinRange(t[1], t[2], 100)
    tEnd = tEnd[2:end] # skip first time step to since propagator won't do anything

    r0_P = x0_P[1:3]
    v0_P = x0_P[4:6]

    dv1_pro = []
    dv1_retro = []
    dv2_pro = []
    dv2_retro = []

    traj_E = []
    traj_prograde = []
    traj_retrograde = []

    for j in 1:length(tEnd)

        # propagate evader orbit to specified time
        x_E_t, x_E_u = propagate_2Body(x₀_E, tEnd[j], mu, 1.0)
        rt_E = x_E_u[end][1:3]
        vt_E = x_E_u[end][4:6]

        # get lambert solution between pursuer state at t to evader state at tStart[end]
        v1_pro, v2_pro = lambertbattin(r0_P, rt_E, mu, "pro", tEnd[j])
        v1_retro, v2_retro = lambertbattin(r0_P, rt_E, mu, "retro", tEnd[j])
        
        # compute dv magnitudes
        push!(dv1_pro, norm(v0_P - v1_pro))
        push!(dv1_retro, norm(v0_P - v1_retro))
        push!(dv2_pro, norm(vt_E - v2_pro))
        push!(dv2_retro, norm(vt_E - v2_retro))

        x_E = mapreduce(permutedims, vcat, x_E_u)
        push!(traj_E, x_E)
        
        # generate output for prograde trajectories        
        x₀_P_lambert = [r0_P; v1_pro]
        x_P_lambert_t, x_P_lambert_u = propagate_2Body(x₀_P_lambert, tEnd[j], mu)
        x_P_lambert = mapreduce(permutedims, vcat, x_P_lambert_u) 
        push!(traj_prograde, x_P_lambert)

        # generate output for retrograde trajectories
        x₀_P_lambert = [r0_P; v1_retro]
        x_P_lambert_t, x_P_lambert_u = propagate_2Body(x₀_P_lambert, tEnd[j], mu)
        x_P_lambert = mapreduce(permutedims, vcat, x_P_lambert_u) 
        push!(traj_retrograde, x_P_lambert)

    end

    generateLambertPlots(x0_P, traj_E, traj_prograde, traj_retrograde, tEnd, dv1_pro, dv1_retro, dv2_pro, dv2_retro, "varyTF")
    return [x0_P, traj_E, traj_prograde, traj_retrograde, tEnd, dv1_pro, dv1_retro, dv2_pro, dv2_retro, "varyTF"]
end 


function generateLambertPlots(fixedTraj, variableTraj, traj_prograde, traj_retrograde, tVals, dv1_pro, dv1_retro, dv2_pro, dv2_retro, source="varyT0")

    dm = ["pro", "retro"]
    text_offset = (0,10)
    for i in 1:length(dm)

        # plot trajectories generated by lambert solutions
        # initialize figure
        fig = Figure() 
        ax = Axis3(fig[1, 1], 
            xlabel = "X (km)", ylabel = "Y (km)", zlabel = "Z (km)") 
        # lines!( fixedTraj[:,1], fixedTraj[:,2], fixedTraj[:,3]; linewidth = 2, label = "Evader" ) 
        # scatter!( fixedTraj[1,1], fixedTraj[1,2], fixedTraj[1,3]; marker = :circle, markersize = 10, color = :black ) 
        # text!( fixedTraj[1,1], fixedTraj[1,2], fixedTraj[1,3]; text = "P IC", color = :gray, offset = text_offset, align = (:center, :bottom) ) 

        for j in 1:length(variableTraj)
            lines!(ax, variableTraj[j][:,1], variableTraj[j][:,2], variableTraj[j][:,3]; linewidth = 2, color = :black)
        end

        if dm[i] == "pro"
            ax.title = "Prograde Lambert Solutions"
            for j in 1:length(traj_prograde)
                lines!(ax, traj_prograde[j][:,1], traj_prograde[j][:,2], traj_prograde[j][:,3]; linewidth = 2)
            end
            save("plots/lambert_"*source*"_prograde.png", fig)
        else
            ax.title = "Retrograde Lambert Solutions"
            for j in 1:length(traj_retrograde)
                lines!(ax, traj_retrograde[j][:,1], traj_retrograde[j][:,2], traj_retrograde[j][:,3]; linewidth = 2)
            end
            save("plots/lambert_"*source*"_retrograde.png", fig)
        end
    end

    # plot both prograde and retrograde trajectories in the same plot 
    fig = Figure() 
    ax = Axis3(fig[1, 1], 
        xlabel = "X (km)", ylabel = "Y (km)", zlabel = "Z (km)") 
    # lines!( fixedTraj[:,1], fixedTraj[:,2], fixedTraj[:,3]; linewidth = 2, label = "Evader" ) 
    # scatter!( fixedTraj[1,1], fixedTraj[1,2], fixedTraj[1,3]; marker = :circle, markersize = 10, color = :black ) 
    # text!( fixedTraj[1,1], fixedTraj[1,2], fixedTraj[1,3]; text = "P IC", color = :gray, offset = text_offset, align = (:center, :bottom) ) 
    for j in 1:length(tVals)
        lines!(ax, traj_prograde[j][:,1], traj_prograde[j][:,2], traj_prograde[j][:,3]; linewidth = 2, color = Makie.wong_colors()[1])
        lines!(ax, traj_prograde[j][:,1], traj_prograde[j][:,2], traj_prograde[j][:,3]; linewidth = 2, color = Makie.wong_colors()[2])
        lines!(ax, variableTraj[j][:,1], variableTraj[j][:,2], variableTraj[j][:,3]; linewidth = 2, color = :black)
    end
    save("plots/lambert_"*source*"_trajectories.png", fig)

    # plot dv1 as a function of t0
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, Float32.(tVals), Float32.(dv1_pro); label="Prograde")
    scatter!(ax, Float32.(tVals), Float32.(dv1_retro); label="Retrograde")
    if source == "varyT0"
        ax.xlabel = "Time of Maneuver 1 (s)"
    else
        ax.xlabel = "Time of Maneuver 2 (s)"
    end    
    ax.ylabel = "|Δv_1| (km/s)"
    ax.title = "Lambert Solution: Intercept (Maneuver 1 only)"
    axislegend(ax, position = :rt)
    save("plots/lambert_"*source*"_intercept.png", fig)

    # plot dv1 as a function of t0
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, Float32.(tVals), Float32.(dv1_pro); color = :red, label="Prograde dv1")
    scatter!(ax, Float32.(tVals), Float32.(dv2_pro); color = :pink, label="Prograde dv2")
    scatter!(ax, Float32.(tVals), Float32.(dv1_retro); color = :blue, label="Retrograde dv1")
    scatter!(ax, Float32.(tVals), Float32.(dv2_retro); color = :teal, label="Retrograde dv2")
    if source == "varyT0"
        ax.xlabel = "Time of Maneuver 1 (s)"
    else
        ax.xlabel = "Time of Maneuver 2 (s)"
    end    
    ax.ylabel = "|Δv| (km/s)"
    ax.title = "Lambert Solution: Rendezvous (Maneuver 1 and 2, Individual)"
    axislegend(ax, position = :rt)
    save("plots/lambert_"*source*"_rendezvous_individual.png", fig)

    # plot sum of dv1 and dv2 as a function of t0
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, Float32.(tVals), Float32.(dv1_pro) + Float32.(dv2_pro); label="Prograde")
    scatter!(ax, Float32.(tVals), Float32.(dv1_retro) + Float32.(dv2_retro); label="Retrograde")
    if source == "varyT0"
        ax.xlabel = "Time of Maneuver 1 (s)"
    else
        ax.xlabel = "Time of Maneuver 2 (s)"
    end    
    ax.ylabel = "|Δv_1| + |Δv_2| (km/s)"
    ax.title = "Lambert Solution: Rendezvous (Maneuver 1 and 2, Total)"
    axislegend(ax, position = :rt)
    save("plots/lambert_"*source*"_rendezvous.png", fig)
end 

export plotTrajAnimation
function plotTrajAnimation(vary_output)
    fixed_traj = vary_output[1]
    variable_traj  = vary_output[2]
    traj_prograde = vary_output[3]
    traj_retrograde = vary_output[4]
    tVals  = vary_output[5]
    dv1_pro = vary_output[6]
    dv1_retro = vary_output[7]
    dv2_pro = vary_output[8]
    dv2_retro = vary_output[9]

    fig = Figure()
    ax = Axis3(fig[1, 1], 
        xlabel = "X (km)", ylabel = "Y (km)", zlabel = "Z (km)")
    N = length(tVals)

    record(fig, "plots/trajAnimation.mp4", 1:1:N) do i
        # plot pursuer and evader trajectories 
        lines!(ax, x_P[:,1], x_P[:,2], x_P[:,3]; linewidth = 2, color = :green)
        lines!(ax, x_E[:,1], x_E[:,2], x_E[:,3]; linewidth = 2, color = :red)

        # plot lambert trajectories
    end


    # record(fig, "plots/trajAnimation.mp4", 1:1:N) do i
    #     _x = LinRange(xmin[i], xmax[i], 200)
    #     _y = LinRange(ymin[i], ymax[i], 200)
    #     hm[1] = _x # update x coordinates
    #     hm[2] = _y # update y coordinates
    #     hm[3] = mandelbrot.(_x, _y') # update data
    #     autolimits!(ax) # update limits
    #     # yield() -> not required with record
    # end
end


export testLambert
function testLambert(r1_t0, r2_t0, mu, t1, t2, dm)
    # propagate orbits to t1 and t2
    r1_t1_t, r1_t1_u = propagate_2Body(r1_t0, t1, mu)
    r1_t1 = r1_t1_u[end]
    
    r2_t0_withDV = applyDV(r2_t0)
    r2_t2_t, r2_t2_u = propagate_2Body(r2_t0_withDV, t2, mu)
    r2_t2 = r2_t2_u[end]

    tof = t2 - t1

    v1, v2 = lambertbattin(r1_t1[1:3], r2_t2[1:3], mu, dm, tof)
    dv1 = norm(r1_t1[4:6] - v1)
    dv2 = norm(r2_t2[4:6] - v2)
    return dv1, dv2
end

export applyDV
function applyDV(r2_t0, dvMag=0/1000)
    r2_t0_new = r2_t0
    
    dir = rand(3)
    dir = dir / norm(dir)
    dv = dvMag * dir

    r2_t0_new[4:6] = r2_t0_new[4:6] + dv

    return r2_t0_new
end

export computeInclinationChange
function computeInclinationChange(v1, Δi)
    kep2 = cart2kep(v1, mu)
    kep2[3] += Δi
    v2 = kep2cart(kep2, mu)[4:6]
    Δi_mag = sqrt(v1^2 + v2^2 - 2*v1*v2*cos(Δi))
    return Δi_mag * (v2 - v1)
end