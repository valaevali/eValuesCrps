library(testthat)

rlang::env_unlock(env = asNamespace('testthat'))
rlang::env_binding_unlock(env = asNamespace('testthat'))
assign('rstudio_tickle', function(){}, envir = asNamespace('testthat'))
rlang::env_binding_lock(env = asNamespace('testthat'))
rlang::env_lock(asNamespace('testthat'))

if (Sys.getenv("NOT_CRAN") == "true") {
  test_check("evalues")
}
