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

function plotCAndApaths( cpath, apath )

    fig = figure("sim2a",figsize=(8,8))
    plot(cpath[:,1],"g",linewidth = 2)
    plot(apath[:,1],"b",linewidth = 2)
    legend("Consumption","Assets");
    xlabel("Age");
    title("Time path of income and consumption")

end
