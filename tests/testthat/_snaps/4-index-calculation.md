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
      ind[, 1:6]
    Output
        year       est       lwr      upr  log_est        se
      1 2003 209276.67 151512.82 289062.8 12.25141 0.1647923
      2 2004 316827.74 244649.37 410300.7 12.66611 0.1319066
      3 2005 355482.83 274277.58 460730.5 12.78123 0.1323169
      4 2007  91541.24  66923.46 125214.7 11.42454 0.1598193
      5 2009 150187.38 109088.68 206769.8 11.91964 0.1631267
      6 2011 274230.39 212808.07 353380.9 12.52172 0.1293789
      7 2013 274541.78 204421.87 368713.9 12.52286 0.1504709
      8 2015 326960.94 241466.64 442725.6 12.69760 0.1546505
      9 2017 151502.42 110227.70 208232.5 11.92836 0.1622751

