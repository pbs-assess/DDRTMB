# Function to plot pairs plots with histogram as used in 2018 P cod assessment
# Originally by Steve Martell, adapted by Chris Grandin and Robyn Forrest

plot.marg <- function(xx,
                      breaks = "sturges",
                      ex.factor = 1.0,
                      specs = 0,
                      ...){
  ## xx - a list(p = samples, p1 = prior param 1, p2 = prior param 2,
  ##  fn = prior distribution)
  ## specs = type of prior
  post.no.plot <- hist(as.matrix(xx$p),
                       breaks = breaks,
                       plot = FALSE)
  xvals <- seq(min(post.no.plot$breaks) / ex.factor,
               max(post.no.plot$breaks) / ex.factor,
               length = 1000)
  pd <- xx$fn(xvals, xx$p1, xx$p2)
  z <- cbind(xvals, pd)

  #set xlim for non-uniform priors
  if(specs==0) {
    xlim <- c(min(xvals), max(xvals))
  }else xlim <- c(0.45*min(xvals), 1.75*max(xvals))

  ss <- hist(as.matrix(xx$p),
             prob = TRUE,
             breaks = breaks,
             main = xx$nm,
             xlab = "",
             cex.axis = 1.2,
             xlim = xlim,
             ylab = "",
             ...)
  func <- function(x){xx$fn(x, xx$p1, xx$p2)}
  ## Plot prior
  curve(func,
        xlim[1],
        xlim[2],
        xlab = "",
        ylab = "",
        col = "black",
        lwd = 2,
        add = TRUE)
  ## Plot MPD
  abline(v = xx$mle,
         lwd = 2,
         lty = 2,
         col = 2)
}
