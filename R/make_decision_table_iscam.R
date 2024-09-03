#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAKE DECISION TABLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make_decision_table_iscam
# Called from run-model-v1.R

# Author: Robyn Forrest
# Date created:  August 2, 2024
# Last Modified: August 3, 2024

# This function takes the outputs of the ISCAM projection model and calculates probabilities
# of various metrics for each TAC
# The only difference between this and make_decision_table is that
# iscam did not output B20xx/B0 so these lines are commented out

# arguments:
# proj_outputs = the list object returned by run_projections()
# npyr = number of projection years

# Returns:
# A decision table

make_decision_table_iscam <- function(proj_output,npyr){

  nyr <-  dat$nyr # actual final historical year
  tac <- pfc$tac.vec

  for(i in seq_along(tac)){
    raw <- proj_output[[i]]

    # Get the probability of each stock status quantity
    tmp1 <-  data.frame(tac[i],
                          mean(raw[,paste0("B",nyr+2,"B",nyr+1)]<1),
                          mean(raw[,paste0("B",nyr+2,"Bavg")]<1),
                          mean(raw[,paste0("B",nyr+2,"Bmin")]<1),
                          mean(raw[,paste0("B",nyr+2,"08Bmsy")]<1),
                          mean(raw[,paste0("B",nyr+2,"04Bmsy")]<1),
                          #mean(raw[,paste0("B",nyr+2,"B0")]<1),
                          mean(raw[,paste0("F",nyr+1,"F",nyr)]>1),
                          mean(raw[,paste0("F",nyr+1,"Favg")]>1),
                          mean(raw[,paste0("F",nyr+1,"Fmsy")]>1))
    # colnames
    cols <- c("TAC",
              paste0("P(B",nyr+2,"<B",nyr+1,")"),
              paste0("P(B",nyr+2,"<Bavg)"),
              paste0("P(B",nyr+2,"<Bmin)"),
              paste0("P(B",nyr+2,"<08Bmsy)"),
              paste0("P(B",nyr+2,"<04Bmsy)"),
              #paste0("P(B",nyr+2,"<B0)"),
              paste0("P(F",nyr+1,">F",nyr,")"),
              paste0("P(F",nyr+1,">Favg)"),
              paste0("P(F",nyr+1,">Fmsy)"))

    if(npyr>1){
        # if more than one projection year, add stock status for subsequent years
        for(ii in 2:npyr){
          tmp2 <- data.frame(mean(raw[,paste0("B",nyr+1+ii,"B",nyr+1)]<1),
                            mean(raw[,paste0("B",nyr+1+ii,"Bavg")]<1),
                            mean(raw[,paste0("B",nyr+1+ii,"Bmin")]<1),
                            mean(raw[,paste0("B",nyr+1+ii,"08Bmsy")]<1),
                            mean(raw[,paste0("B",nyr+1+ii,"04Bmsy")]<1),
                            #mean(raw[,paste0("B",nyr+1+ii,"B0")]<1),
                            mean(raw[,paste0("F",nyr+ii,"F",nyr)]>1),
                            mean(raw[,paste0("F",nyr+ii,"Favg")]>1),
                            mean(raw[,paste0("F",nyr+ii,"Fmsy")]>1))
          # colnames
          cols <- c(cols,
                    paste0("P(B",nyr+1+ii,"<B",nyr+1,")"),
                    paste0("P(B",nyr+1+ii,"<Bavg)"),
                    paste0("P(B",nyr+1+ii,"<Bmin)"),
                    paste0("P(B",nyr+1+ii,"<08Bmsy)"),
                    paste0("P(B",nyr+1+ii,"<04Bmsy)"),
                    #paste0("P(B",nyr+1+ii,"<B0)"),
                    paste0("P(F",nyr+ii,">F",nyr,")"),
                    paste0("P(F",nyr+ii,">Favg)"),
                    paste0("P(F",nyr+ii,">Fmsy)"))
        tmp1 <- cbind(tmp1,tmp2)
      } # end for ii
    }# end if

   if(i==1){
     dtable <- tmp1
   }else{
     dtable <- rbind(dtable,tmp1)
   }
  }# end for i (tac loop)

  colnames(dtable) <- cols
  return(as.data.frame(dtable))

} # end function

