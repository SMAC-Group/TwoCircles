load("Simulation study/simulation3b.RData")

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
                               "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
                               "\\usepackage{amssymb}"))
## I need the amssymb package because I use \mathcal and \mathbb
tikz("Simulation study/tex/simu3b.tex", width = 7, height = 3.5, standAlone = TRUE,
     packages = c("\\usepackage{tikz}",
                  "\\usepackage[active,tightpage,psfixbb]{preview}",
                  "\\PreviewEnvironment{pgfpicture}",
                  "\\setlength\\PreviewBorder{0pt}",
                  "\\usepackage{amssymb}"))


# Graphics of functions
par(mfrow = c(1,2),  mar=c(4,3.5,2,2), mgp = c(3, 1, 0))
mu.max = 20
mu.length = 15
delta = seq(from = 0, to = mu.max, length.out = mu.length)

taille.axes = 1
taille.leg = 1

# Define colors
n.fun = 20
coleur = topo.colors(n.fun)
index = c(1,3,6,8,11,13)
coleur = c(coleur[index[1]], "deepskyblue2", "forestgreen", "chocolate")
index = c(1,2,3,4)


plot(NA, ylim = c(0,1), xlim = range(delta), main = "dg", xlab = "x", ylab = "y" , cex.axis = 0.7, cex = 0.7, cex.lab = 0.8, ann = FALSE, axes = FALSE)
grid(col = "grey85")
mtext("Rejection frequency of H$_0$", side = 2, line = 2, cex=taille.leg, col = "grey30")
mtext("$\\delta$", side = 1, line = 2, cex=taille.leg, col = "grey30")
box(col = "black")
title("\\textbf{(a)}")
axis(2, col = "black", labels = TRUE, outer=FALSE , cex.axis = taille.axes)
axis(1,  col = "black", labels = TRUE, outer=FALSE, cex.axis = taille.axes)

abline(h = 0.05, lwd = 2)

k = 1
lines(delta, out$res_rao[k,]/100, col = coleur[2], lty = 1)
lines(delta, out$res_dixon[k,]/100, col = coleur[3], lty = 2)
lines(delta, out$res_WW[k,]/100, col = coleur[4], lty = 4)

plot(NA, ylim = c(0,1), xlim = range(delta), main = "dg", xlab = "x", ylab = "y" , cex.axis = 0.7, cex = 0.7, cex.lab = 0.8, ann = FALSE, axes = FALSE)
grid(col = "grey85")
mtext("Rejection frequency of H$_0$", side = 2, line = 2, cex=taille.leg, col = "grey30")
mtext("$\\delta$", side = 1, line = 2, cex=taille.leg, col = "grey30")
box(col = "black")
title("\\textbf{(b)}")
axis(2, col = "black", labels = TRUE, outer=FALSE , cex.axis = taille.axes)
axis(1,  col = "black", labels = TRUE, outer=FALSE, cex.axis = taille.axes)

abline(h = 0.05, lwd = 2)

k = 2
lines(delta, out$res_rao[k,]/100, col = coleur[2], lty = 1)
lines(delta, out$res_dixon[k,]/100, col = coleur[3], lty = 2)
lines(delta, out$res_WW[k,]/100, col = coleur[4], lty = 4)

legend("topleft",c("Rao","Dixon","Wheeler-Watson"),
       col = coleur[c(2,3,4)], lty = c(1,2,4), bty = "n", bg = "white", cex = 1)


dev.off()

