source("renv/activate.R")

# https://stackoverflow.com/questions/4996090/how-to-disable-save-workspace-image-prompt-in-r
utils::assignInNamespace(
  "q", 
  function(save = "no", status = 0, runLast = TRUE) 
  {
    .Internal(quit(save, status, runLast))
  }, 
  "base"
)
