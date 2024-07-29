# This code modified from the 2018 Pacific Cod Res Doc
# https://github.com/pbs-assess/pacific-cod-2018/blob/master/R/figures-mcmc-diagnostics.R

make.priors.posts.plot <- function(mc,
                                   ctl,
                                   priors.only = TRUE){
  ## Make a plot of the priors used in the given model overlaid on the
  ##  posterior
  ##
  ## priors.only - Only plot the priors and not posteriors

  f.names <- c(dunif, dnorm, dlnorm, dbeta, dgamma)
  post.names <- names(mc)

  prior.specs <- as.data.frame(ctl$params)
  ## Remove fixed parameters
  prior.specs <- prior.specs[prior.specs$phz > 0,]
  ## Remove upper and lower bound, and phase information, but keep initial
  ##  value
  prior.specs <- prior.specs[, -c(2:4)]
  prior.names <- rownames(prior.specs)

  ## Add the q parameters to the prior specs table
  ## NOTE that RTMB reports them all in natural space
  ## even though prior is in log space
  q.params <- ctl$surv.q
  num.q.params <- ncol(q.params)
  q.specs <- lapply(1:num.q.params,
                    function(x){
                      c(q.params[2, x],
                        q.params[1, x],
                        q.params[2, x],
                        q.params[3, x])
                    })
  q.specs <- as.data.frame(do.call(rbind, q.specs))
  rownames(q.specs) <- paste0("log_q", 1:num.q.params)
  colnames(q.specs) <- colnames(prior.specs)
  prior.specs <- rbind(prior.specs, q.specs)

  ## Remove log part for q's with uniform priors
  ## Remove log part for q's with uniform priors
  j <- prior.specs
  non.q <- j[-grep("q", rownames(j)),]
  non.q <- non.q %>%
    rownames_to_column() %>%
    as.tibble()
  q <- j[grep("q", rownames(j)),]
  q <- q %>%
    rownames_to_column() %>%
    mutate(rowname = if_else(!prior,
                             gsub("log_", "", rowname),
                             rowname))
  prior.specs <- as.data.frame(rbind(non.q, q))
  rownames(prior.specs) <- prior.specs$rowname
  prior.specs <- prior.specs %>% select(-rowname)
  rownames(prior.specs)[rownames(prior.specs) == "steepness"] = "h"

  # RTMB: For q with lognormal prior
  #  Need to convert the q posterior to log space to match prior
  for(i in 1:num.q.params){
    if(q$prior[i]==1){
      mc[,paste0("q",i)] <- log(mc[,paste0("q",i)])
    }
  }

  # If priors only, need to change dimensions of plot to reduce white space
  # Also need to remove priors with uniform dist
  if(priors.only==TRUE){
    prior.specs <-  prior.specs %>%  dplyr::filter(prior > 0)
    post.names <- row.names(prior.specs)
  }

  n.side <- get.rows.cols(length(post.names))
  par(mfrow = n.side,
      oma = c(2, 3, 1, 1),
      mai = c(0.3, 0.4, 0.3, 0.2))

  for(i in 1:length(post.names)){
    specs <- prior.specs[i,]
    prior.fn <- f.names[[as.numeric(specs[2] + 1)]]
    xx <- list(p = mc[,i],
               p1 = as.numeric(specs[3]),
               p2 = as.numeric(specs[4]),
               fn = prior.fn,
               nm = rownames(prior.specs)[i])

    xx$nm <- get.latex.name(xx$nm)

    if(priors.only){
      func <- function(x){xx$fn(x, xx$p1, xx$p2)}
      if(specs[2] == 0){
        ## Uniform, plot from p1-1 to p2+1
        ## curve(func,
        ##       from = xx$p1 - 1,
        ##       to = xx$p2 + 1,
        ##       xlab = "",
        ##       ylab = "",
        ##       col = "black",
        ##       lwd = 2)
      }else if(specs[2] == 1){
        ## Normal, plot from -(p1-p2*4) to (p1+p2*4)
        curve(func,
              from = xx$p1 - 4 * xx$p2,
              to = xx$p2 + 2 * xx$p2,
              xlab = "",
              ylab = "",
              col = "black",
              lwd = 2)
        title(xx$nm)
      }else{
        curve(func,
              xlab = "",
              ylab = "",
              col = "black",
              lwd = 2)
        title(xx$nm)
      }
    }else{
      plot.marg(xx,
                breaks = "sturges",
                col = "wheat",
                specs = prior.specs[i,2])
    }
  }
}

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

make.traces.plot <- function(model,
                             axis.lab.freq = 200){
  ## Make trace plots for all paramaters from the mcmc output
  ## axis.lab.freq - the frequency of x-axis labelling

  if(class(model) == model.lst.class){
    model <- model[[1]]
    if(class(model) != model.class){
      stop("The structure of the model list is incorrect.")
    }
  }

  mc <- model$mcmc$params.est
  ## Remove some of them
  ## mc <- mc[, -grep("ro", colnames(mc))]
  mc <- mc[, -grep("rinit", colnames(mc))]
  mc <- mc[, -grep("rbar", colnames(mc))]
  mc <- mc[, -grep("bo", colnames(mc))]
  mc <- mc[, -grep("msy", colnames(mc))]
  mc <- mc[, -grep("ssb", colnames(mc))]

  n.side <- get.rows.cols(ncol(mc))
  par(mfrow = n.side,
      oma = c(2, 3, 1, 1),
      mai = c(0.2, 0.4, 0.3, 0.2))

  for(param in 1:ncol(mc)){
    mcmc.trace <- as.matrix(mc[,param])
    name <- colnames(mc)[param]
    name <- get.latex.name(name)
    plot(mcmc.trace,
         main = name,
         type = "l",
         ylab = "",
         xlab = "",
         axes = FALSE)
    box()
    at <- labels <- seq(0,
                        nrow(mc),
                        axis.lab.freq)
    axis(1,
         at = at,
         labels = labels)
    axis(2)
  }
}

make.autocor.plot <- function(model,
                              ylim = c(-1,1)){
  ## Plot the autocorrelation of estimated parameters

  if(class(model) == model.lst.class){
    model <- model[[1]]
    if(class(model) != model.class){
      stop("The structure of the model list is incorrect.")
    }
  }

  mc <- model$mcmc$params.est
  ## Remove some of them
  ## mc <- mc[, -grep("ro", colnames(mc))]
  mc <- mc[, -grep("rinit", colnames(mc))]
  mc <- mc[, -grep("rbar", colnames(mc))]
  mc <- mc[, -grep("bo", colnames(mc))]
  mc <- mc[, -grep("msy", colnames(mc))]
  mc <- mc[, -grep("ssb", colnames(mc))]

  n.side <- get.rows.cols(ncol(mc))
  par(mfrow = n.side,
      oma = c(2, 3, 1, 1),
      mai = c(0.2, 0.4, 0.3, 0.2))

  for(param in 1:ncol(mc)){
    mcmc.autocor <- as.matrix(mc[,param])
    name <- colnames(mc)[param]
    name <- get.latex.name(name)
    autocorr.plot(mcmc.autocor,
                  lag.max = 100,
                  main = name,
                  auto.layout = FALSE,
                  ylim = ylim)
  }
}

autocorr.plot <- function(x,
                          lag.max,
                          auto.layout = TRUE,
                          ask,
                          ylim = c(-1,1),
                          ...){
  ## autocorr.plot from coda package, but the package source had the
  ##  ylab = "Autocorrelation" for all plots and no way to override it.
  ## That caused latex to place the plot in landscape mode which was ugly.
  if (missing(ask)){
    ask <- if(is.R()) {
      dev.interactive()
    }else{
      interactive()
    }
  }
  oldpar <- NULL
  on.exit(par(oldpar))
  if(auto.layout)
    oldpar <- par(mfrow = set.mfrow(Nchains = nchain(x),
                                    Nparms = nvar(x)))
  if(!is.mcmc.list(x))
    x <- mcmc.list(as.mcmc(x))
  for(i in 1:nchain(x)) {
    xacf <- if(missing(lag.max))
              acf(as.ts.mcmc(x[[i]]),
                  plot = FALSE)
            else
              acf(as.ts.mcmc(x[[i]]),
                  lag.max = lag.max,
                  plot = FALSE)
    for(j in 1:nvar(x)){
      plot(xacf$lag[-1, j, j],
           xacf$acf[-1, j, j],
           type = "h",
           ylab = "",
           xlab = "Lag",
           ylim = ylim, ...)
      title(paste0(varnames(x)[j],
                  ifelse(is.null(chanames(x)),
                         "",
                         ":"),
                  chanames(x)[i]))
      if(i == 1 & j == 1)
        oldpar <- c(oldpar, par(ask = ask))
    }
  }
  invisible(x)
}

as.ts.mcmc <- function(x, ...){
  ## as.ts.mcmc was copied from coda package source, to fulfill
  ##  autocorr.plot requirement

  x <- as.mcmc(x)
  y <- ts(x,
          start = start(x),
          end = end(x),
          deltat = thin(x))
  attr(y, "mcpar") <- NULL
  y
}

make.pairs.plot <- function(mc){

  ## Plot the pairs scatterplots for the estimated parameters
  panel.hist <- function(x, ...){
    usr    <- par("usr");
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h      <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y      <- h$counts; y <- y/max(y)
    rect(breaks[-nB],
         0,
         breaks[-1],
         y,
         col="wheat",
         cex=0.75,
         ...)
  }

  mc <- as.data.frame(mc)
  ## Remove some of them
  ## mc <- mc[, -grep("ro", colnames(mc))]
  # mc <- mc[, -grep("rinit", colnames(mc))]
  # mc <- mc[, -grep("rbar", colnames(mc))]
  # mc <- mc[, -grep("bo", colnames(mc))]
  # mc <- mc[, -grep("msy", colnames(mc))]
  # mc <- mc[, -grep("ssb", colnames(mc))]

  c.names <- colnames(mc)
  latex.names <- NULL
  for(param in 1:ncol(mc)){
    name <- c.names[param]
    latex.names <- c(latex.names, get.latex.name(name))
  }

  our.panel.smooth <- function (x, y, col = "#00000060", bg = NA, pch = par("pch"),
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)
  {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok))
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
        col = col.smooth, ...)

    r <- round(cor(x, y), digits=2)
    txt <- sprintf("%.2f", r)
    par(usr = c(0, 1, 0, 1))
    graphics::text(0.98, 0.9, txt, pos = 2, col = "grey5")
  }

  pairs(mc,
        labels = latex.names,
        cex.labels = 1.0,
        pch = ".",
        upper.panel = our.panel.smooth,
        diag.panel = panel.hist,
        lower.panel = our.panel.smooth,
        gap = 0.0)
}

get.rows.cols <- function(num){
  ## Returns a vector of length 2 representing the number of rows and columns
  ##  to use to pack a plot in a grid.
  if(num <= 64 && num > 49){
    if(num <= 56){
      nside <- c(8,7)
    }else{
      nside <- c(8,8)
    }
  }else if(num <= 49 && num > 36){
    if(num <= 42){
      nside <- c(7,6)
    }else{
      nside <- c(7,7)
    }
  }else if(num <= 36 && num > 25){
    if(num <= 30){
      nside <- c(6,5)
    }else{
      nside <- c(6,6)
    }
  }else if(num <= 25 && num > 16){
    if(num <= 20){
      nside <- c(5,4)
    }else{
      nside <- c(5,5)
    }
  }else if(num <= 16 && num > 9){
    if(num <= 12){
      nside <- c(4,3)
    }else{
      nside <- c(4,4)
    }
  }else if(num <=  9 && num > 4){
    if(num <= 6){
      nside <- c(3,2)
    }else{
      nside <- c(3,3)
    }
  }else if(num <=  4 && num > 1){
    if(num == 2){
      nside <- c(2,1)
    }else{
      nside <- c(2,2)
    }
  }else{
    nside <- c(1,1)
  }
  return(nside)
} # end get.rows.cols

get.latex.name <- function(name, addToQ = 0){
  ## Return a pretty version of the parameter name found in variable 'name'
  ##
  ## addToQ - an integer to the parameter name for the q's. This is necessary
  ##  because iscam sets the q parameter names to 1, 2, 3... regardless of the
  ##  gear number. i.e. if gear 1 is a trawl fishery and gear 2 is a survey,
  ##  iscam will call q1 the survey gear. We must add 1 to it to get q2 to
  ##  accurately portray the survey gear number
  if(name == "ro") return(expression("R"[0]))
  if(name == "rbar") return(expression(bar("R")))
  if(name == "rinit") return(expression(bar("R")[init]))
  if(name == "m") return(expression("M"))
  if(name == "bo") return(expression("B"[0]))
  if(name == "vartheta") return(expression(vartheta))
  if(name == "rho") return(expression(rho))
  if(name == "bmsy") return(expression("B"[MSY]))
  if(name == "msy") return(expression("MSY"))
  if(name == "fmsy") return(expression("F"[MSY]))
  if(name == "umsy") return(expression("U"[MSY]))
  if(name == "ssb") return(expression("SSB"))
  if(name == "sel1") return(expression(hat(a)[1]))
  if(name == "selsd1") return(expression(hat(gamma)[1]))
  if(name == "sel2") return(expression(hat(a)[2]))
  if(name == "selsd2") return(expression(hat(gamma)[2]))
  if(name == "sel3") return(expression(hat(a)[3]))
  if(name == "selsd3") return(expression(hat(gamma)[3]))
  if(name == "sel4") return(expression(hat(a)[4]))
  if(name == "selsd4") return(expression(hat(gamma)[4]))
  if(name == "sel5") return(expression(hat(a)[5]))
  if(name == "selsd5") return(expression(hat(gamma)[5]))
  if(name == "log_ro") return(expression("ln(R"[0]*")"))
  if(name == "h") return(expression("h"))
  if(name == "m1") return(expression("M"[1]))
  if(name == "m2") return(expression("M"[2]))
  if(name == "log_m") return(expression("ln(M)"))
  if(name == "log_rbar") return(expression("ln("*bar("R")*")"))
  if(name == "log_rinit") return(expression("ln("*bar("R")[init]*")"))

  if(length(grep("^q[1-9]+$", name))){
    digit <- as.numeric(sub("^q([1-9]+)$", "\\1", name))
    return(substitute("q"[digit], list(digit = digit)))
  }

  if(length(grep("^log_q[1-9]+$", name))){
    digit <- as.numeric(sub("^log_q([1-9]+)$", "\\1", name))
    return(substitute("ln(q"[digit]*")", list(digit = digit)))
  }

  NULL
} # end get.latex.name
