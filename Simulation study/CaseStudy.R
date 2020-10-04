library(circular)

control1 = c(75, 80,95,130,170,210)/360*2*pi
control2 = c(75,80)/360*2*pi
control3 = 80/360*2*pi
experimental1 = c(10,50,55,65,90,285,325,355)/360*2*pi
experimental2 = c(55,285)/360*2*pi


#control1 = c(74, 78,80,82,95,130,170,210)/360*2*pi
#control1 = c(76, 78,80,82,95,130,170,210)
#experimental = c(10,50,55,55,65,90,285,285,325,355)

control.circular = circular(control, type = c("angles"),
units = c("degrees"),
template = c("none"),
modulo = c("asis"),
zero = 0,
rotation = c("counter"),
name = "control")

experimental.circular = circular(experimental, type = c("angles"),
units = c("degrees"),
template = c("none"),
modulo = c("asis"),
zero = 0,
rotation = c("counter"),
name = "Experimental")

empty.circular = circular(NA, type = c("angles"),
units = c("degrees"),
template = c("none"),
modulo = c("asis"),
zero = 0,
rotation = c("counter"),
name = "Experimental")




# Load tikz -> lib
library(tikzDevice)

# Latex options
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
    "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
    "\\usepackage{amssymb}","\\usepackage{amsfonts}"))

width = 5
height = 5

# Title of latex file + latex packages needed
tikz("Simulation study/tex/case_study.tex", width = width, height = height, standAlone = TRUE,
    packages = c("\\usepackage{tikz}",
                "\\usepackage[active,tightpage,psfixbb]{preview}",
                 "\\PreviewEnvironment{pgfpicture}",
                 "\\setlength\\PreviewBorder{0pt}",
                 "\\usepackage{amssymb}",
				"\\usepackage{bm}","\\usepackage{amsthm}","\\usepackage{amsbsy}"
				,"\\usepackage{amsbsy}"
				,"\\usepackage{amsbsy}"
				,"\\usepackage{amsfonts}"))

par(mfrow = c(1,1))
plot(empty.circular, cex = 1.3, ylim = c(-1, 1.2), xlim = c(-1, 1))


points(cos(control1), sin(control1), col=1, pch = 16 , cex = 1.8)
points(cos(experimental1), sin(experimental1), col=1, pch = 21, cex = 1.8, bg = "white")

delta = 0.1
points((1 + delta)*cos(control2), (1 + delta)*sin(control2), col=1, pch = 16 , cex = 1.8)
points((1 + delta)*cos(experimental2), (1 + delta)*sin(experimental2), col=1, pch = 21, cex = 1.8, bg = "white")

points((1 + 2*delta)*cos(control3), (1 + 2*delta)*sin(control3), col=1, pch = 16 , cex = 1.8)

legend("topleft",c("Control","Experimental"), pch = c(16,21), pt.cex = c(2,2),
cex = 1, bty = "n", bg = "white", box.col = "white")


dev.off()
