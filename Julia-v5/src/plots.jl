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
            size	= (1200,800)
            )

        plot!(cpath[:,i],
            lc = "green",
            linewidth = 2,
            label = "Consumption"
            )

        plot!(apath[:,i],
            lc = "blue",
            linewidth = 2,
            label = "Assets"
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
