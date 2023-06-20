
#RC_PR checks and returns the recall and precision of an estimator

RC_PR <- function(coef_vector,
                 active_preds){
  rc = sum(which(coef_vector != 0) %in% active_preds)/length(active_preds)
  if(sum(coef_vector != 0))
    pr = sum(which(coef_vector != 0) %in% active_preds)/sum(coef_vector != 0) else 
      pr = 0
  
  return(
    list(rc = rc, pr =pr)
  )
  }

