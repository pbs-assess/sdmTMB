# tests of arguments in sdmTMBcontrol()

SEED <- 123
set.seed(SEED)

test_that("nlminb improves model fit with sdmTMBcontrol()", {
  local_edition(2)
  skip_if_not_installed("INLA")
  d <- subset(pcod, year >= 2009)
  pcod_spde <- make_mesh(d, c("X", "Y"), cutoff = 15)
  m_1 <- sdmTMB(density ~ 1,
                 data = d, time = "year", mesh = pcod_spde,
                 control = sdmTMBcontrol(nlminb_loops = 1),
                 family = tweedie(link = "log"), spatiotemporal = "RW"
  )
  m_2 <- sdmTMB(density ~ 1,
                 data = d, time = "year", mesh = pcod_spde,
                 control = sdmTMBcontrol(nlminb_loops = 2),
                 family = tweedie(link = "log"), spatiotemporal = "RW"
  )
  expect_true(max(m_1$gradients) > max(m_2$gradients))
})
