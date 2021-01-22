######### Table Manipulation Function ########

pvalFormat <- function(p.values, method = 'none', replace = FALSE){
  ## Formats p-values for reports, can report adjusted pvalues
  ##    Inputs:
  ##       - p.value: numeric p-value
  ##       - method: pvalue adjustment, passed to p.adjust.methods
  ##       - replace: if TRUE, replaces p-values with their adjusted value
  ##    Outputs:
  ##       - out: formatted p-value

  p.values <- suppressWarnings(as.numeric(p.values))
  out      <- rep(NA, length(p.values))
  sig      <- p.adjust(p.values, method)
  if(replace) p.values <- sig

  for(i in 1:length(p.values)){
    if(is.na(p.values[i])){out[i] <- NA}else{
      if(p.values[i] >= .001){
        out[i] <- paste('$', formatC(p.values[i], format = 'f', digits = 3), '$', sep = '')
      }

      if(p.values[i] < .001){
        out[i] <- '< $0.001$'
      }

      if(sig[i] > 0.01 & sig[i] <= 0.05){
        out[i] <- paste(out[i], '*', sep = '')
      }

      if(sig[i] > 0.001 & sig[i] <= 0.01) {
        out[i] <- paste(out[i], '**', sep = '')
      }

      if(sig[i] <= 0.001){
        out[i] <- paste(out[i], '***', sep = '')
      }}
  }

  out[is.na(out)] <- '-'
  return(out)
}


