.onAttach <- function(...) {
  packageStartupMessage("drtmle: TMLE with doubly robust inference")
  packageStartupMessage(
    "Version: ",
    utils::packageDescription("drtmle")$Version
  )
}
