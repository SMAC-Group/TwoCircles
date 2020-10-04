n = 10
m = 9
x = get_critical_values(n = n, m = m, test = "rao")
tab = table(x$dist)

p = c(0, cumsum(as.numeric(tab/sum(tab))))

cols = c("black", "grey60") #hcl(h = seq(15, 375, length = 3), l = 65, c = 100)[1:2]

# Graphics of functions


options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
                               "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
                               "\\usepackage{amssymb}"))
## I need the amssymb package because I use \mathcal and \mathbb
tikz("Simulation study/tex/fig1.tex", width = 7, height = 3.5, standAlone = TRUE,
     packages = c("\\usepackage{tikz}",
                  "\\usepackage[active,tightpage,psfixbb]{preview}",
                  "\\PreviewEnvironment{pgfpicture}",
                  "\\setlength\\PreviewBorder{0pt}",
                  "\\usepackage{amssymb}"))
par(mfrow = c(1,2),  mar=c(3,4,2,2), mgp = c(1, 1, 0))
main =  "Rao"
plot(NA, xlim=c(-2, 17), ylim= c(0, 1), lwd=2, ylab=" ",
     xlab = " ", main = main, cex.main = 1.2)
grid()
mtext("$P(x)$", side = 2, line = 3)
mtext("$T_1$", side = 1, line = 2)
xx = seq(from = -2, to = 20, length.out = 10^4)
yy = pnorm(xx, m, sqrt(m/2))
lines(xx, yy, col = cols[1], lwd = 1.5, lty = 2)

xx = c(-10,as.numeric(names(tab)),100)
for (i in 1:length(p)){
  lines(c(xx[i], xx[i+1]), c(p[i], p[i]), col = cols[2])
  if (i > 1 && i < length(p)){
    lines(c(xx[i+1],xx[i+1]), c(p[i], p[i+1]), lty = 1, col = cols[2])
  }
  #points(xx[i+1],p[i],pch=21,cex=0.75,col=cols[2], bg = "white")
  #points(xx[i],p[i],pch=16,cex=0.75,col=cols[2])
}


x = get_critical_values(n = n, m = m, test = "dixon")
tab = table(x$dist)

p = c(0, cumsum(as.numeric(tab/sum(tab))))

legend("topleft", c("Asymptotic", "Exact"), lwd = 1.5, lty = 2:1, col = cols, bty = "n")
#cols = hcl(h = seq(15, 375, length = 3), l = 65, c = 100)[1:2]


main =  "Dixon"
plot(NA, xlim=c(0, 70), ylim= c(0, 1), lwd=2, ylab=" ",
     xlab = " ", main = main, cex.main = 1.2)
grid()
mtext("$P(x)$", side = 2, line = 3)
mtext("$T_2$", side = 1, line = 2)
xx = seq(from = -2, to = 100, length.out = 10^4)
yy = pnorm(xx, 3*m, 4*sqrt(m))
lines(xx, yy, col = cols[1], lwd = 1.5, lty = 2)

xx = c(-10,as.numeric(names(tab)),200)
for (i in 1:length(p)){
  lines(c(xx[i], xx[i+1]), c(p[i], p[i]), col = cols[2])
  if (i > 1 && i < length(p)){
    lines(c(xx[i+1],xx[i+1]), c(p[i], p[i+1]), lty = 1, col = cols[2])
  }
  #points(xx[i+1],p[i],pch=21,cex=0.75,col=cols[2], bg = "white")
  #points(xx[i],p[i],pch=16,cex=0.75,col=cols[2])
}
dev.off()

