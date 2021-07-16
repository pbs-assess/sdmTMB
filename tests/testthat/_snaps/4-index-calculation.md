# get_index(), get_index_sims(), and get_cog() work

    Code
      m
    Output
      Spatial model fit by ML ['sdmTMB']
      Formula: density ~ 0 + as.factor(year)
      Time column: "year"
      SPDE: pcod_spde
      Data: pcod
      Family: tweedie(link = 'log')
                          coef.est coef.se
      as.factor(year)2003     2.76    0.37
      as.factor(year)2004     3.09    0.35
      as.factor(year)2005     3.09    0.36
      as.factor(year)2007     1.95    0.37
      as.factor(year)2009     2.27    0.37
      as.factor(year)2011     2.88    0.36
      as.factor(year)2013     3.00    0.35
      as.factor(year)2015     3.03    0.36
      as.factor(year)2017     2.35    0.37
      
      Matern range parameter: 3.82
      Dispersion parameter: 13.69
      Spatial SD: 11.38
      Spatiotemporal SD: 7.52
      ML criterion at convergence: 6526.833
      
      See ?tidy.sdmTMB to extract these values as a data frame.

---

    Code
      ind
    Output
        year       est       lwr      upr  log_est        se max_gradient bad_eig
      1 2003 209277.00 151513.00 289063.4 12.25141 0.1647925 4.085145e-05   FALSE
      2 2004 316826.47 244648.36 410299.1 12.66611 0.1319067 4.085145e-05   FALSE
      3 2005 355481.95 274276.86 460729.4 12.78123 0.1323169 4.085145e-05   FALSE
      4 2007  91540.99  66923.25 125214.4 11.42454 0.1598195 4.085145e-05   FALSE
      5 2009 150187.86 109089.00 206770.6 11.91964 0.1631269 4.085145e-05   FALSE
      6 2011 274230.96 212808.47 353381.7 12.52173 0.1293790 4.085145e-05   FALSE
      7 2013 274541.85 204421.91 368714.0 12.52286 0.1504710 4.085145e-05   FALSE
      8 2015 326960.55 241466.30 442725.1 12.69759 0.1546506 4.085145e-05   FALSE
      9 2017 151502.37 110227.66 208232.4 11.92836 0.1622751 4.085145e-05   FALSE

