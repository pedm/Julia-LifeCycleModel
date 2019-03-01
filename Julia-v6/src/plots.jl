# To learn more about plotting:
# https://gist.github.com/gizmaa/7214002


function plotCpath( cpath )
    plot(cpath,
        linewidth = 2,
        title="Time path of consumption",
        xlabel="Age",
        size	= (1200,800)
        )
    gui()
end

function plotApath( apath, borrowCon )

    plot(apath,
        linewidth = 2,
        title = "Time path of assets",
        xlabel="Age",
        size = (1200, 800)
        )

    plot!(borrowCon,
        lc = "red",
        linewidth = 5,
        label = "borrowing con"
        )

    gui()

end

function plotYAndCpaths( ypath, cpath )

    for i = 1:2

        plot(ypath[:,i],
            lc = "red",
            linewidth = 2,
            xlabel = "Age",
            label = "Income",
            title = string("Time path of income and consumption - individual ", i),
            size	= (1200,800)
            )

        plot!(cpath[:,i],
            lc = "green",
            linewidth = 2,
            label = "Consumption"
            )

        gui()
    end

    # fig = figure("sim2",figsize=(8,8))
    # plot(ypath[:,2],"r",linewidth = 2)
    # plot(cpath[:,2],"g",linewidth = 2)
    # legend("Income","Consumption");
    # xlabel("Age");
    # title("Time path of income and consumption - individual 2")

end

function plotYCAndApaths( ypath, cpath, apath )

    for i = 1:2

        plot(ypath[:,i],
            lc = "red",
            linewidth = 2,
            xlabel = "Age",
            label = "Income",
            title = string("Time path of individual ", i),
            size	= (1200,800),
            hover = ypath[:,i]
            )

        plot!(cpath[:,i],
            lc = "green",
            linewidth = 2,
            label = "Consumption",
            hover = cpath[:,i]
            )

        plot!(apath[:,i],
            lc = "blue",
            linewidth = 2,
            label = "Assets",
            hover = apath[:,i]
            )
        gui()
    end

end

# Subplots function is similar to that in matlab:

# # fig, axes = subplots(2, 2, figsize=(16,6))
# fig = figure("test",figsize=(8,8))
# subplot(221) # Create the 1st axis of a 2x2 arrax of axes
# plot(apath)
# title("Assets")
#
# subplot(222)
# plot( collect(1:T), cpath)
# # plot(ypath)
# # legend("C", "Y")
# title("Consumption")
#
# subplot(223)
# plot(collect(1:T), ypath)
# title("Income")

function plot_V(t = 5, ixY = numPointsY)
    EVt = EV[t, :, ixY]
    plot(Agrid[t, :], EVt, label = "EV", marker = :circle, size = (1200, 800), title = "EV and V at time t = $t and ixY = $ixY", hover = EVt)
    for ixtYtr = 1:numPointsYTrans
        Vt = V[t, :, ixY, ixtYtr]
        plot!(Agrid[t, :], Vt, label = "ixtYtr = $ixtYtr", marker = :circle, hover = Vt)
    end
    gui()
end

function plot_policyA1(t = 5, ixY = numPointsY)
    plot(size = (1200, 800), title = "Policy Function for Assets at time t = $t and ixY = $ixY")
    for ixtYtr = 1:numPointsYTrans
        A1 = policyA1[t, :, ixY, ixtYtr]
        plot!(Agrid[t, :], A1, label = "ixtYtr = $ixtYtr", marker = :circle, hover = A1)
    end
    gui()
end

function plot_EV_over_time()
    # Recall EV dimensions are T+1, numPointsA, numPointsY
    t= 5
    plot(Agrid[t, :], EV[t, :, numPointsY], label = "t = $t", marker = :circle, size = (1200, 800), title = "EV")
    for t = 15:10:55
        plot!(Agrid[t, :], EV[t, :, numPointsY], label = "t = $t", marker = :circle)
    end
    gui()
end
