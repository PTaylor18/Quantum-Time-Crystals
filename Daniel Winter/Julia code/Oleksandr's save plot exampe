using Plots

y(x) = x .^2 # apply function to each element in the list x
x = collect(0:0.1:10)

plt = plot(x, y(x), line = (:solid, 2), xaxis = ("x"), yaxis = ("y"), label = ("y(x)"), legend = true);
display(plt)

macro Name(arg)
   string(arg)
end

function save_plot(name, plot, label)
    cd(dirname(@__FILE__));
    dir = pwd();
    println("Saving plots to the directory in $dir","\\plots")
    strlabel = string(label);
    mkpath(string(dir,"\\plots"))
    pngfilename = string(dir,"\\plots\\",name,"_",strlabel,".png")
    pdffilename = string(dir,"\\plots\\",name,"_",strlabel,".pdf")
    #epsfilename = string(dir,"\\plots\\",name,"_",strlabel,".eps")
    savefig(plot, pngfilename)
    savefig(plot, pdffilename)
    #savefig(plot, epsfilename)
end

label="parabola";
plot_name = @Name plot
save_plot(plot_name, plt, label)


println("Finished")
