function plotCpath( cpath )
    fig = figure("cpath",figsize=(8,8))
    plot(cpath, linewidth = 2)
    xlabel("Age");
    title("Time path of consumption")
end

function plotApath( apath, borrowCon )
    fig = figure("apath",figsize=(8,8))
    plot(apath, linewidth = 2)
    plot(borrowCon, "--r", linewidth = 2)
    # legend("1", "2", "Borrowing con");
    xlabel("Age")
    title("Time path of assets")

end

function plotCAndApaths( cpath, apath )

    display(Plots.plot([1:length(cpath[:, 1])], cpath[:,1],linewidth = 2))

end
