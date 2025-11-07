# Get TMB parameter list

Get TMB parameter list

## Usage

``` r
get_pars(object)
```

## Arguments

- object:

  Fit from
  [`sdmTMB()`](https://pbs-assess.github.io/sdmTMB/reference/sdmTMB.md)

## Value

A named list of parameter values

## Examples

``` r
fit <- sdmTMB(present ~ 1, data = pcod_2011, family = binomial(), spatial = "off")
pars <- get_pars(fit)
names(pars)
#>  [1] "ln_H_input"          "b_j"                 "b_j2"               
#>  [4] "bs"                  "ln_tau_O"            "ln_tau_Z"           
#>  [7] "ln_tau_E"            "ln_kappa"            "thetaf"             
#> [10] "gengamma_Q"          "logit_p_extreme"     "log_ratio_mix"      
#> [13] "ln_phi"              "ln_tau_V"            "rho_time_unscaled"  
#> [16] "ar1_phi"             "re_cov_pars"         "re_b_pars"          
#> [19] "b_rw_t"              "omega_s"             "zeta_s"             
#> [22] "epsilon_st"          "b_threshold"         "b_epsilon"          
#> [25] "ln_epsilon_re_sigma" "epsilon_re"          "b_smooth"           
#> [28] "ln_smooth_sigma"    
```
