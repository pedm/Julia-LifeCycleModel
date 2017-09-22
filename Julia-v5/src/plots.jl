# To learn more about plotting:
# https://gist.github.com/gizmaa/7214002

function plotCpath( cpath )
    fig = figure("cpath",figsize=(8,8))
    plot(cpath, linewidth = 2)
    xlabel("Age");
    title("Time path of consumption")
end

function plotApath( apath, borrowCon )
    fig = figure("apath",figsize=(8,8))
    plot(apath, lineWidth = 2)
    plot(borrowCon, "--r", linewidth = 2)
    # legend("1", "2", "Borrowing con");
    xlabel("Age")
    title("Time path of assets")

end

function plotYAndCpaths( ypath, cpath )

    fig = figure("sim1",figsize=(8,8))
    plot(ypath[:,1],"r",linewidth = 2)
    plot(cpath[:,1],"g",linewidth = 2)
    legend("Income","Consumption");
    xlabel("Age");
    title("Time path of income and consumption - individual 1")

    fig = figure("sim2",figsize=(8,8))
    plot(ypath[:,2],"r",linewidth = 2)
    plot(cpath[:,2],"g",linewidth = 2)
    legend("Income","Consumption");
    xlabel("Age");
    title("Time path of income and consumption - individual 2")

end

function plotYCAndApaths( ypath, cpath, apath )

    fig = figure("sim2a",figsize=(8,8))
    plot(ypath[:,1],"r",linewidth = 2)
    plot(cpath[:,1],"g",linewidth = 2)
    plot(apath[:,1],"b",linewidth = 2)
    xlabel("Age");
    title("Time path of income and consumption - individual 1")


    fig = figure("sim2b",figsize=(8,8))
    plot(ypath[:,2],"r",linewidth = 2)
    plot(cpath[:,2],"g",linewidth = 2)
    plot(apath[:,2],"b",linewidth = 2)
    xlabel("Age");
    title("Time path of income and consumption - individual 2")

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
