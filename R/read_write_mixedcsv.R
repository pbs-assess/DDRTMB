#functions originally written by Michael Folkes for fun. 



#this still needs to be adapted to double nesting lists
#readMixedcsv <- function(filename){
#
#  data.tmp <- readLines(filename)
#  
#  data.tmp <- trimws(data.tmp)
#  df <- data.frame(comment.lines = grep("^#", data.tmp))
#  
#  df$data.rows.length <- diff(c(df$comment.lines, length(data.tmp)+1))-1
#  
#  comments <- data.tmp[df$comment.lines]
#  comments <- gsub("^# *", "", comments)
#  
#  data.list <- apply(df, 1, function(index){
#    index <- as.list(index)
#    if(index$data.rows.length>0){
#    read.table(text=data.tmp, skip = index$comment.lines, nrows = index$data.rows.length, sep = ",")
#    }
#  })
#  names(data.list) <- comments
#  data.list
#}#END readMixedcsv


writeMixedcsv <- function(x, filename="test.csv"){

  file.create(filename)

  invisible(
  lapply(1:length(x), function(index){
    commentline <- paste("#", names(x)[index])
    write.table(commentline, file = filename, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    if(is.list(x[[index]])&!is.data.frame(x[[index]])){
      xx<-x[[index]]
       lapply(1:length(xx), function(subindex){
         commentline <- paste("#", names(xx)[subindex])
         write.table(commentline, file = filename, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
         if(is.data.frame(xx[[subindex]])){
            write.table(xx[subindex], file = filename, append = TRUE, col.names = TRUE, row.names = FALSE, sep=",", quote = FALSE)
          }else{
            write.table(xx[subindex], file = filename, append = TRUE, col.names = FALSE, row.names = FALSE, sep=",", quote = FALSE)
          }
         
        })
    }else{
      if(is.data.frame(x[[index]])){
        write.table(x[index], file = filename, append = TRUE, col.names = TRUE, row.names = FALSE, sep=",", quote = FALSE)
      }else{
        write.table(x[index], file = filename, append = TRUE, col.names = FALSE, row.names = FALSE, sep=",", quote = FALSE, qmethod="double")
      }
      
    }
    
  })
 )

}#END writeMixedcsv


