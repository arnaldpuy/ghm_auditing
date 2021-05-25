## ----setup, include=FALSE--------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dev = "tikz")


## ----initial, results="hide", message=FALSE, warning=FALSE-----------------------

# PRELIMINARY STEPS ---------------------------------------------------------------

# Install and load packages in one go
loadPackages <- function(x) {
  for (i in x) {
    if (!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}

# Load the packages
loadPackages(c(
  "data.table", "tidyverse", "sensobol", "wesanderson",
  "triangle", "scales", "cowplot", "fitdistrplus",
  "parallel", "doParallel", "foreach", "sp", "sf", "Rfast",
  "raster", "rworldmap", "countrycode"))

# Custom theme
theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent",
                                           color = NA),
          legend.key = element_rect(fill = "transparent",
                                    color = NA),
          legend.position = "top",
          strip.background = element_rect(fill = "white"))
}

# set checkpoint
dir.create(".checkpoint")
library("checkpoint")

checkpoint("2021-02-22", 
           R.version ="4.0.3", 
           checkpointLocation = getwd())


## ----irrigated_areas, cache=TRUE-------------------------------------------------

# PLOT UNCERTAINTY IN IRRIGATED AREAS -----------------------------------------------

# Read data compiled by Meier
irrigated.areas <- fread("meier.csv")

# Arrange and drop Oceania and 0's
dt <- melt(irrigated.areas, measure.vars = 4:9) %>%
  .[!Continent == "Oceania"] %>%
  .[!value == 0]

# List to plot
continent_list <- list(c("Africa", "Americas"), c("Asia", "Europe"))

# Plot
gg <- list()
for (i in 1:length(continent_list)) {
  gg[[i]] <- ggplot(dt[Continent %in% continent_list[[i]]], 
                    aes(reorder(Country, value), value)) +
    geom_point(stat = "identity", aes(color = variable)) +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    coord_flip() +
    scale_color_manual(
      name = "Dataset",
      values = c(
        "yellowgreen", "seagreen4", "magenta3",
        "sienna3", "turquoise2", "khaki3"
      )
    ) +
    labs(
      x = "",
      y = "Irrigated area (ha)"
    ) +
    facet_wrap(~Continent, scales = "free_y") +
    theme_AP()
}


## ----plot_irrigated_areas, cache=TRUE, dependson="irrigated_areas", fig.height=6.5, fig.width=5.7, warning=FALSE, message=FALSE, fig.cap="Irrigated area estimates produced for Africa and the Americas by the FAO-GMIA [@Siebert2013], the IWMI--GIAM [@Thenkabail2009], the GRIPC [@Salmon2015], Meier's map [@Meier2018], Aquastat [@FAO2016b] and FAOSTAT [@FAO2017a]. The data has been retrieved from @Meier2018"----

gg[[1]]


## ----plot_irrigated_areas2, cache=TRUE, dependson="irrigated_areas", fig.height=6.5, fig.width=5.7, warning=FALSE, message=FALSE, fig.cap="Irrigated area estimates produced for Asia and Europe by the FAO-GMIA [@Siebert2013], the IWMI--GIAM [@Thenkabail2009], the GRIPC [@Salmon2015], Meier's map [@Meier2018], Aquastat [@FAO2016b] and FAOSTAT [@FAO2017a]. The data has been retrieved from @Meier2018"----

gg[[2]]


## ----plot_kc, cache=TRUE, fig.height=2.5, fig.width=5.5, fig.cap="The crop coefficient. a) Evolution of $k_c$ values for salt cedar. The solid black line is the mean $k_c$, the red dots show the individual measured values (adapted from Figure 10 in @Nichols2004), and the vertical, dashed black line marks the values selected for the global sensitivity analysis presented in Figure~1c and d. b) Oasis effect on $k_c$ values (adapted from Figure 46 in @Allen1998)."----

# K_C VALUES FOR SALT CEDAR ----------------------------------------------------

# Read in dataset
kc_evolution_cedar <- fread("kc_evolution_cedar.csv")
kc_evolution_point <- fread("kc_evolution_point.csv")

# Retrieve minimum and maximmum values for month 5
cedar.min.max <- kc_evolution_point[x > 5 & x < 6]

# Plot
a <- ggplot(kc_evolution_cedar, aes(x = x, y = y)) +
  geom_line() +
  geom_point(
    data = kc_evolution_point, aes(x = x, y = y),
    color = "red", alpha = 0.4
  ) +
  scale_color_discrete(name = "") +
  annotate("text", x = 5.78, y = 0.8, label = "$k_c$ dev") +
  annotate("text", x = 8, y = 1.4, label = "$k_c$ med") +
  annotate("text", x = 9, y = 0.9, label = "$k_c$ late") +
  geom_vline(xintercept = 5, lty = 2) +
  theme_AP() +
  labs(x = "Month", y = "$k_c$")

# Read in data to show oasis effect
oasis <- fread("oasis_effect.csv")

# Plot and merge
b <- ggplot(oasis, aes(x = x, y = y, group = name, linetype = name)) +
  geom_line() +
  scale_linetype_discrete(name = "") +
  theme_AP() +
  labs(x = "Width of irrigation area (m)", y = "") +
  theme(legend.position = c(0.6, 0.5))

cowplot::plot_grid(a, b, ncol = 2, labels = "auto", rel_widths = c(0.45, 0.55))


## ----kc_only, cache=TRUE, dependson="plot_kc", fig.height=2.2, fig.width=2.4-----

# PLOT K_C VALUES FOR SALT CEDAR ONLY ------------------------------------------

a


## ----calculations, cache=TRUE, dependson="plot_kc"-------------------------------

# DEFINE SETTINGS --------------------------------------------------------------

N <- 2^12
R <- 10^3

# LISTED FUNCTIONS -------------------------------

et0_fun <- list(

  "pt" = function(alpha, delta, gamma, A, k_c = 1)
    k_c * (alpha * ((delta * A) / (gamma + delta))),

  "pm" = function(delta, A, gamma, T_a, w, v, k_c = 1)
    k_c * ((0.408 * delta * A + gamma * (900 / (T_a + 273)) * w * v) /
                    delta + gamma * (1 + 0.34 * w))
)

mat_transform <- function(mat, method, et_c = FALSE) {

  if (method == "pt") {
    mat[, "alpha"] <- qunif(mat[, "alpha"], 1.26 + 1.26 * - 0.1, 1.26 + 1.26 * 0.1)
    mat[, "delta"] <- qunif(mat[, "delta"], 0.21 + 0.21 * - 0.005, 0.21 + 0.21 * 0.005)
    mat[, "gamma"] <- qunif(mat[, "gamma"], 0.059 + 0.059 * - 0.001, 0.059 + 0.059 * 0.001)
    mat[, "A"] <- qunif(mat[, "A"], 350 + 350 * - 0.15, 350 + 350 * 0.15)

  } else if (method == "pm") {
    mat[, "delta"] <- qunif(mat[, "delta"], 0.21 + 0.21 * - 0.005, 0.21 + 0.21 * 0.005)
    mat[, "gamma"] <- qunif(mat[, "gamma"], 0.059 + 0.059 * - 0.001, 0.059 + 0.059 * 0.001)
    mat[, "A"] <- qunif(mat[, "A"], 350 + 350 * - 0.15, 350 + 350 * 0.15)
    mat[, "T_a"] <- qunif(mat[, "T_a"], 27 + 27 * - 0.01, 27 + 27 * 0.01)
    mat[, "w"] <- qunif(mat[, "w"], 2.5 + 2.5 * - 0.05, 2.5 + 2.5 * 0.05)
    mat[, "v"] <- qunif(mat[, "v"], 2.49 + 2.49 * - 0.04, 2.49 + 2.49 * 0.04)
  }

  if (et_c == TRUE) {
    mat[, "k_c"] <- qunif(mat[, "k_c"], min(cedar.min.max$y), max(cedar.min.max$y))
  }
  return(mat)
}

params.shared <- c("delta", "gamma", "A")
params.pt <- c(params.shared, "alpha")
params.pm <- c(params.shared, "T_a", "v", "w")

# COMPUTATION OF ET_0 ------------------------------------------------------------

y <- list()

for (i in list("pt", "pm")) {

  if(i == "pt") {

    mat <- sobol_matrices(N = N, params = params.pt)
    mat <- mat_transform(mat = mat, method = i, et_c = FALSE)
    y[[i]] <- et0_fun[[i]](alpha = mat[, "alpha"],
                      delta = mat[, "delta"],
                      gamma = mat[, "gamma"],
                      A = mat[, "A"])

  } else if (i == "pm") {

    mat <- sobol_matrices(N = N, params = params.pm)
    mat <- mat_transform(mat = mat, method = i, et_c = FALSE)
    y[[i]] <- et0_fun[[i]](delta = mat[, "delta"],
                      A = mat[, "A"],
                      gamma = mat[, "gamma"],
                      T_a = mat[, "T_a"],
                      w = mat[, "w"],
                      v = mat[, "v"])
  }
}

# COMPUTATION OF ET_C ------------------------------------------------------------

y.kc <- list()

for(i in list("pt", "pm")) {

  if(i == "pt") {

    mat <- sobol_matrices(N = N, params = c(params.pt, "k_c"))
    mat <- mat_transform(mat = mat, method = i, et_c = TRUE)
    y.kc[[i]] <- et0_fun[[i]](alpha = mat[, "alpha"],
                           delta = mat[, "delta"],
                           gamma = mat[, "gamma"],
                           A = mat[, "A"],
                           k_c = mat[, "k_c"])

  } else if (i == "pm") {

    mat <- sobol_matrices(N = N, params = c(params.pm, "k_c"))
    mat <- mat_transform(mat = mat, method = i, et_c = TRUE)
    y.kc[[i]] <- et0_fun[[i]](delta = mat[, "delta"],
                           A = mat[, "A"],
                           gamma = mat[, "gamma"],
                           T_a = mat[, "T_a"],
                           w = mat[, "w"],
                           v = mat[, "v"],
                           k_c = mat[, "k_c"])
  }
}

# FUNCTION TO EXTRACT UNCERTAINTY ------------------------------------------------

extract_unc <- function(out, N = N) {
  value_output <- unlist(lapply(out, function(x) x[1:N]), use.names = FALSE)
  da <- data.table(value_output)[, method:= rep(names(et0_fun), each = N)] %>%
    .[, method:= toupper(method)]
  gg <- ggplot(da, aes(value_output, fill = method)) +
    geom_histogram(alpha = 0.3, color = "black", position = "identity") +
    scale_fill_manual(name = "$ET_0$ formula",
                      values = wes_palette(n = 2, name = "Chevalier1")) +
    theme_AP() +
    theme(legend.position = "none")
  return(gg)
}

# PLOT UNCERTAINTY ---------------------------------------------------------------

plots.unc <- lapply(list(y, y.kc), function(x) extract_unc(out = x, N = N))

plots.unc[[1]] <- plots.unc[[1]] + labs(x = "$ET_0$ (mm d$^{-1}$)", y = "Counts") +
  scale_x_continuous(limits = c(0, 800))
plots.unc[[2]] <- plots.unc[[2]] + labs(x = "$ET_c$ (mm d$^{-1}$)", y = "Counts") +
  scale_x_continuous(limits = c(0, 800))

etc <- plots.unc[[2]] + theme(legend.position = "top") +
  scale_fill_manual(labels = c("FAO-56 Penman Monteith", "Priestley-Taylor"), 
                    values = wes_palette(n = 2, name = "Chevalier1"))
  # For later

# SOBOL' INDICES -----------------------------------------------------------------

params.pt.plot <- c("$\\Delta$", "$\\gamma$", "$A$", "$\\alpha$")
params.pm.plot <- c("$\\Delta$", "$\\gamma$", "$A$", "$T_a$", "$v$", "$w$")
first <- "jansen"

# ET0 --------------------------

ind <- list()

for (i in list("pt", "pm")) {

  if (i == "pt") {
    params <- params.pt.plot

  } else if (i == "pm") {
    params <- params.pm.plot
  }

  ind[[i]] <- sobol_indices(Y = y[[i]], N = N, params = params,
                            first = first, boot = TRUE, R = R)
  ind[[i]]$results[, method:= i]
}

# ET0 --------------------------

ind.kc <- list()

for( i in list("pt", "pm")) {

  if (i == "pt") {
    params <- c(params.pt.plot, "$k_c$")

  } else if (i == "pm") {
    params <- c(params.pm.plot, "$k_c$")
  }

  ind.kc[[i]] <- sobol_indices(Y = y.kc[[i]], N = N, params = params,
                               first = first, R = R, boot = TRUE)
  ind.kc[[i]]$results[, method:= i]
}

# FUNCTION TO EXTRACT SOBOL' INDICES ---------------------------------------------

extract_sobol <- function(data) {
  out <- lapply(data, function(x) x$results) %>%
    rbindlist(.) %>%
    .[, method:= toupper(method)] %>%
    .[, parameters:= factor(parameters, levels = c("$\\Delta$", "$\\gamma$",
                                                   "$A$", "$T_a$", "$w$", "$v$",
                                                   "$\\alpha$", "$k_c$"))] %>%
    ggplot(., aes(parameters, original, fill = sensitivity)) + 
    geom_bar(stat = "identity", 
             position = position_dodge(0.6), color = "black") + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
    labs(x = "", y = "Sobol' index") + 
    ggplot2::scale_fill_discrete(name = "Sobol' indices", 
                                 labels = c(expression(S[italic(i)]), expression(T[italic(i)]))) + 
    theme_AP() +
    theme(legend.position = "none") +
    facet_grid(~ method,
               scales = "free_x",
               space = "free_x")
  return(out)
}

# PLOT SOBOL  --------------------------------------------------------------------

plots.ind <- lapply(list(ind, ind.kc), function(x) extract_sobol(data = x))


## ----plot, cache=TRUE, dependson="calculations", fig.height=4.8, fig.width=5, message=FALSE, warning=FALSE, fig.cap="Monte-Carlo based, global uncertainty and sensitivity analysis of the Penman-Monteith (PM) and the Priestley-Taylor (PT) methods. The red bar indicates the first-order effect $S_i$, i.e. the proportion of variance conveyed by each model input. The blue bar reflects the total-order effect $T_i$, which includes the first-order effect of the parameter plus the effect derived of its interaction with the rest. The parameters were described with the probability distributions shown in Table~1. See the Supplementary Materials for a technical explanation of global uncertainty and sensitivity analysis. a) Empirical distribution of $ET_0$. b) Sobol' indices. c) Empirical distribution of $ET_c$. We described $k_c$ as $U(0.47, 1.86)$ to reproduce the range reported at month 5 by @Nichols2004 for developing salt cedars in the Bosque del Apache, New Mexico (see Figure~3a). d) Sobol' indices."----

# MERGE ALL PLOTS ----------------------------------------------------------------

# Arrange legend
legend_ua <- cowplot::get_legend(plots.unc[[1]] + theme(legend.position = "top"))
legend_sa <- cowplot::get_legend(plots.ind[[1]] + theme(legend.position = "top"))
all_legend <- plot_grid(legend_ua, legend_sa, ncol = 2, rel_heights = c(0.01, 0.01))

# Arrange plots
top <- plot_grid(plots.unc[[1]], plots.ind[[1]], ncol = 2, labels = c("a", "b"))
bottom <- plot_grid(plots.unc[[2]], plots.ind[[2]], ncol = 2, labels = c("c", "d"))
all_plot <- plot_grid(top, bottom, ncol = 1)

# Merge
plot_grid(all_legend, all_plot, ncol = 1, rel_heights = c(0.1, 0.9))


## ----efficiencies_usa_africa, cache=TRUE, fig.height=5, fig.width=4.2, message=FALSE, warning=FALSE, fig.cap="Data on irrigation efficiency. a) Distribution of project efficiencies ($E_p$) in USA according to @Solley1998, calculated as Total water withdrawal / consumptive water use. The vertical  dashed line is the $E_p$ value used by @Doll2002 to characterize the irrigation efficiency of USA. b) Distribution of project efficiencies ($E_p$) for Africa reported by @FAO1997. c) Distribution of field application efficiencies ($E_a$) in surface and sprinkler irrigation systems according to @Bos1990. ``Surface`` includes furrow, basin and border irrigation systems. The dashed vertical lines mark the $E_a$ point estimates selected by @Rohwer2007. d) Distribution of conveyance efficiencies ($E_c$) in surface and pressurized (sprinkler) irrigation systems according to @Bos1990. The vertical, dashed black lines show the $E_c$ point estimates used by @Rohwer2007."----

# PLOT EFFICIENCIES --------------------------------------------------------------

# USA data
usa.dt <- fread("usa_efficiency.csv")
usa.dt <- usa.dt[, Efficiency:= consumptive.use / total.withdrawal]

a <- ggplot(usa.dt, aes(Efficiency)) +
  geom_histogram() +
  geom_vline(xintercept = 0.6, lty = 2) +
  labs(x = "$E_p$", y = "Counts") +
  theme_AP()

# FAO 1997 (Irrigation potential in Africa)
fao_dt <- fread("fao_1997.csv")
fao_dt <- fao_dt[, Efficiency:= Efficiency / 100]

b <- ggplot(fao_dt, aes(Efficiency)) +
  geom_histogram() +
  labs(x = "$E_p$", y = "") +
  theme_AP()

# Bos and Nugteren data
bos.dt <- fread("bos.dt.csv")
col_names <- colnames(bos.dt)[2:7]
setnames(bos.dt, col_names, paste("$", col_names, "$", sep = ""))

# Create data set with E_a values as defined by Rohwer
dt_e_a <- data.table("Type" = c("Sprinkler", "Surface"),
                     "Value" = c(0.75, 0.6))

c <- ggplot(bos.dt, aes(`$E_a$`)) +
  geom_histogram(bins = 15) +
  geom_vline(data = dt_e_a, aes(xintercept = Value), lty = 2) +
  facet_grid(~ Type) +
  labs(x = "$E_a$", y = "Counts") +
  theme_AP()

# Create data set with E_c values as defined by Rohwer
dt_e_c <- data.table("Type" = c("Surface", "Pressurized"),
                     "Value" = c(0.7, 0.95))

bos.dt.copy <- copy(bos.dt)

d <- bos.dt.copy[, Type:= ifelse(Type == "Surface", "Surface", "Pressurized")] %>%
  ggplot(., aes(`$E_c$`)) +
  geom_histogram(bins = 15) +
  geom_vline(data = dt_e_c, aes(xintercept = Value), lty = 2) +
  facet_grid(~ Type) +
  labs(x = "$E_c$", y = "Counts") +
  theme_AP()

# Merge all plots
top <- cowplot::plot_grid(a, b, ncol = 2, labels = "auto")
bottom <- cowplot::plot_grid(c, d, ncol = 1, labels = c("c", "d"))
cowplot::plot_grid(top, bottom, ncol = 1, rel_heights = c(0.35, 0.65))


## ----efficiencies_scale, cache=TRUE, dependson="efficiencies_usa_africa", fig.height=3, fig.width=4.2, fig.cap = "Partial efficiencies as a function of scale according to the data compiled by @Bos1990. $E_d$ stands for distribution efficiency and refers to the water delivery from the source to a specific outlet. $E_d$ is replaced by @Rohwer2007 with $M_f$."----

# EFFICIENCIES AS A FUNCTION OF SCALE ----------------------------------------------

dt.tmp <- bos.dt[, Scale := ifelse(Irrigated_area < 10000,
  "$<10.000$ ha", "$>10.000$ ha"
)] %>%
  na.omit()

melt(dt.tmp, measure.vars = c("$E_c$", "$E_a$", "$E_d$")) %>%
  ggplot(., aes(value)) +
  geom_histogram() +
  labs(x = "Value", y = "Counts") +
  facet_grid(Scale ~ variable) +
  theme_AP()


## ----settings_factorial, cache=TRUE----------------------------------------------

# DEFINE SETTINGS ------------------------------------------------------------------

N <- 2^14
params <- c("E_a", "E_c", "M_f")
R <- 10^3


## ----matrix_factorial, cache=TRUE, dependson="settings_factorial"----------------

# DEFINE SAMPLE MATRIX -------------------------------------------------------------

mat <- sobol_matrices(N = N, params = params, type = "R")

# Truncated Weibull for E_a
shape <- 3.502469
scale <- 0.5444373
minimum <- 0.14
maximum <- 0.87
Fa.weibull <- pweibull(minimum, shape = shape, scale = scale)
Fb.weibull <- pweibull(maximum, shape = shape, scale = scale)
mat[, "E_a"] <- qunif(mat[, "E_a"], Fa.weibull, Fb.weibull)
mat[, "E_a"] <- qweibull(mat[, "E_a"], shape, scale)

# Truncated Beta for E_c
shape1 <- 5.759496
shape2 <- 1.403552
minimum <- 0.26
maximum <- 0.98
Fa.beta <- pbeta(minimum, shape1 = shape1, shape2 = shape2)
Fb.beta <- pbeta(maximum, shape1 = shape1, shape2 = shape2)
mat[, "E_c"] <- qunif(mat[, "E_c"], Fa.beta, Fb.beta)
mat[, "E_c"] <- qbeta(mat[, "E_c"], shape1, shape2)

# Truncated Weibull for M_f
shape <- 6.844793
scale <- 0.8481904
minimum <- 0.5
maximum <- 0.97
Fa.weibull <- pweibull(minimum, shape = shape, scale = scale)
Fb.weibull <- pweibull(maximum, shape = shape, scale = scale)
mat[, "M_f"] <- qunif(mat[, "M_f"], Fa.weibull, Fb.weibull)
mat[, "M_f"] <- qweibull(mat[, "M_f"], shape, scale)


## ----distributions_factorial, cache=TRUE, dependson="matrix_factorial", fig.height=3,fig.width=4, fig.cap="Empirical and modeled distributions for the uncertainty and sensitivity analysis of Equation 5."----

# PLOT DISTRIBUTIONS ---------------------------------------------------------------

dt <- data.table(mat)
params_plot <- paste("$", params, "$", sep = "")
dt <- setnames(dt, params, params_plot) %>%
  .[, Data:= "Modeled"]

# Filter empirical datasets
dt.ea <- bos.dt[Type == "Surface", `$E_a$`]
dt.ec <- bos.dt[Type == "Surface", `$E_c$`]
dt.mf <- bos.dt[Irrigated_area < 10000, `$E_d$`]

n <- max(length(dt.ea), length(dt.ec), length(dt.mf))

length(dt.ea) <- n
length(dt.ec) <- n
length(dt.mf) <- n

dt.empirical <- cbind(dt.ea, dt.ec, dt.mf) %>%
  data.table() %>%
  .[, Data:= "Empirical"] %>%
  setnames(., colnames(.), colnames(dt))

rbind(dt.empirical, dt[1:N]) %>%
  melt(., measure.vars = params_plot) %>%
  ggplot(., aes(value)) +
  geom_histogram() +
  labs(x = "Value", y = "Counts") +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  facet_grid(Data ~ variable, scales = "free_y") +
  theme_AP()


## ----check_corr, cache=TRUE, dependson=c("distributions_factorial"), fig.height=2, fig.width=4, fig.cap="Scatterplots of $E_a,E_c,M_f$. The topmost and bottomost label facets refer to the $x$ and the $y$ axis respectively."----

# PLOT SCATTERPLOTS TO CHECK FOR CORRELATION BETWEEN E_A, E_C, M_F ---------------

out <- t(combn(colnames(dt.empirical)[1:3], 2))
da <- list()

for (i in 1:nrow(out)) {
  cols <- out[i, ]
  da[[i]] <- cbind(dt.empirical[, ..cols], cols[1], cols[2])
  data.table::setnames(da[[i]], colnames(da[[i]]), c("xvar",
                                                     "yvar", "x", "y"))
}

output <- data.table::rbindlist(da)

ggplot(output, aes(xvar, yvar)) +
  geom_point() +
  facet_wrap(x ~ y) +
  labs(x = "", y = "") +
  theme_AP()


## ----model_factorial, cache=TRUE, dependson="matrix_factorial", fig.height=4, fig.width=4----

# DEFINE MODEL ---------------------------------------------------------------------

y <- mat[, "E_a"] * mat[, "M_f"] * mat[, "E_c"]

# Some statistics
data.table(y)[, .(min = min(y), max = max(y))]

quantile(y, probs = c(0.025, 0.975))


## ----uncertainty_factorial, cache=TRUE, dependson="model_factorial"--------------

# PLOT UNCERTAINTY -----------------------------------------------------------------

a <- plot_uncertainty(Y = y, N = N) +
  labs(x = "Irrigation efficiency", y = "Counts") +
  geom_vline(xintercept = 0.38, lty = 2, color = "red")
a

ep <- a # For later


## ----sensitivity_factorial, cache=TRUE, dependson="model_factorial", fig.height=2.8, fig.width=2.8----

# SENSITIVITY ANALYSIS -------------------------------------------------------------

ind <- sobol_indices(
  Y = y, N = N, params = params_plot, R = R, boot = TRUE,
  first = "jansen"
)
b <- plot_sobol(ind) +
  theme(legend.position = "none")
b


## ----full_sensitivity_factorial, cache=TRUE, dependson=c("model_factorial", "sensitivity_factorial"), fig.height=2.3, fig.width=5.5----

# MERGE PLOTS ----------------------------------------------------------------------

legend <- get_legend(b + theme(legend.position = "top"))
bottom <- cowplot::plot_grid(a, b, ncol = 2, labels = "auto")
cowplot::plot_grid(legend, bottom, ncol = 1, rel_heights = c(0.15, 0.85))


## ----oat_exploration, cache=TRUE, fig.height=2, fig.width=2.2--------------------

# ASSESS THE FRACTION OF THE UNCERTAIN SPACE EXAMINED BY OAT -----------------------

# Formula: Ratio of the hypercube to the hypersphere
oat_exploration <- function(k) pi^((k) / 2) * (0.5)^(k) / gamma(1 + k / 2)

# Check from 1 to 20 dimensions
out <- sapply(1:20, function(x) oat_exploration(x))

# Plot
data.table(k = 1:20, x = out) %>%
  ggplot(., aes(k, x)) +
  geom_point() +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  geom_hline(yintercept = 0.05, lty = 2, color = "red") +
  geom_line() +
  labs(
    y = "Fraction explored",
    x = "$k$"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    legend.key = element_rect(
      fill = "transparent",
      color = NA
    )
  )


## ----functions_maps, cache=TRUE--------------------------------------------------

# FUNCTIONS TO CONVERT LON LAT TO USA STATES AND COUNTRIES -------------------------

# Lon Lat to USA states
# (extracted from https://stackoverflow.com/questions/
# 8751497/latitude-longitude-coordinates-to-state-code-in-r)
lonlat_to_state <- function(pointsDF,
                            states = spData::us_states,
                            name_col = "NAME") {
  ## Convert points data.frame to an sf POINTS object
  pts <- st_as_sf(pointsDF, coords = 1:2, crs = 4326)

  ## Transform spatial data to some planar coordinate system
  ## (e.g. Web Mercator) as required for geometric operations
  states <- st_transform(states, crs = 3857)
  pts <- st_transform(pts, crs = 3857)

  ## Find names of state (if any) intersected by each point
  state_names <- states[[name_col]]
  ii <- as.integer(st_intersects(pts, states))
  state_names[ii]
}

# Lon Lat to countries
coords2country <- function(points) {
  countriesSP <- rworldmap::getMap(resolution = "low")
  pointsSP <- sp::SpatialPoints(points, proj4string = CRS(proj4string(countriesSP)))
  indices <- sp::over(pointsSP, countriesSP)
  indices$ADMIN
}


## ----read_maps, cache=TRUE, dependson="functions_maps"---------------------------

# READ IN RASTERS ------------------------------------------------------------------

# Define parallel computing
n_cores <- floor(detectCores() * 0.75)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Vector with the name of the files
c("fao_gmia.asc", "meier_map.tif", "iwmi_giam.tif")


vec_rasters <- c("fao_gmia.asc", "GRIPC_irrigated_area.asc")

# Load rasters and transform to csv in parallel
out <- foreach(
  i = 1:length(vec_rasters),
  .packages = c(
    "raster", "data.table", "sf",
    "rworldmap", "sp"
  )
) %dopar% {
  rs <- raster(vec_rasters[i])
  dt <- data.table(rasterToPoints(rs,
    fun = function(r) {
      r > 0
    }
  ))
  states_vector <- lonlat_to_state(dt[, 1:2])
  dt[, states := cbind(states_vector)]
  dt[, country := coords2country(dt[, c("x", "y")])]
  setnames(dt, 3, "area")
}

# Stop parallel cluster
stopCluster(cl)


## ----modify_rasters, cache=TRUE, dependson="read_maps"---------------------------

# ARRANGE DATASET -----------------------------------------------------------------

# Read the meier map
meier.map <- fread("meier.map.csv")

# Arrange dataset
names(out) <- c("FAO-GMIA", "GRIPC")
rasters.dt <- rbindlist(out, idcol = "Map")

# Rbind GRIPC, FAO-GMIA and Meier map
all.rasters <- rbind(rasters.dt, meier.map)

# Check differences in global irrigated areas
all.rasters[, sum(area) / 10^6, Map] # Million ha

# Export
fwrite(all.rasters, "all.rasters.csv")

# Coordinates of Uvalde, Texas
# 29.209684, -99.786171.

# keep for later: [x< -99.6 & x > -99.8 & y > 29 & y < 29.5]

# meier map: x = -99.70416, y = 29.34583
rasters.merge <- copy(rasters.dt)
rasters.merge <- rasters.merge[, c("x", "y") := round(.SD, 4), .SDcols = c("x", "y")]
rasters.uvalde <- rbind(
  rasters.merge[x == -99.7083 & y == 29.4583],
  meier.map[x >= -99.71 & x <= -99.70 &
    y >= 29.2 & y <= 29.46]
)


## ----diff_grid, cache=TRUE, dependson="modify_rasters", fig.height=3.5, fig.width=3.5----

# CHECK DIFFERENCES AT THE GRID CELL LEVEL PER CONTINENT ---------------------------

rasters.dt <- rasters.dt[, c("x", "y"):= round(.SD, 4), .SDcols = c("x", "y")]
tmp <- merge(rasters.dt[Map == "FAO-GMIA"], rasters.dt[Map == "GRIPC"], by = c("x", "y"))

# Compute absolute difference at the grid level
tmp <- tmp[, abs:= abs(area.x - area.y)]

# Get continent
tmp <- tmp[, continent:= countrycode(tmp[, country.x],
                                     origin = "country.name",
                                     destination = "continent")]

# Plot
tmp[!continent == "Oceania"] %>%
ggplot(., aes(abs)) +
  geom_histogram() +
  labs(x = "Ha", y = "Counts") +
  facet_wrap(~continent) +
  theme_AP()


## ----kc_wheat, cache=TRUE--------------------------------------------------------

# CHECK THE UNCERTAINTY IN KC COEFFICIENTS FOR WHEAT IN TEXAS -----------------

kc_wheat <- fread("kc_wheat_new.csv")

# Define the time frame of the data
min.day <- 50
max.day <- 53

# Retrieve the filtered data
kc_wheat.dt <- kc_wheat[x > min.day & x < max.day][y > 0]

# Uniform distribution
v.final <- runif(10^4, min = min(kc_wheat.dt$y), max = max(kc_wheat.dt$y))


## ----plot_kc_wheat, cache=TRUE, dependson="kc_wheat", fig.height=2, fig.width=5.5----

# PLOT DATA, EMPIRICAL DISTRIBUTION AND MODELED DISTRIBUTION -----------------------

a <- ggplot(kc_wheat, aes(x = x, y = y)) +
  annotate("rect",
    xmin = min.day, xmax = max.day,
    ymin = 0, ymax = Inf, fill = "red"
  ) +
  geom_point(size = 0.6) +
  geom_smooth() +
  labs(x = "Days after planting", y = "$k_c$") +
  theme_AP()

b <- ggplot(kc_wheat.dt, aes(y)) +
  geom_histogram() +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  labs(x = "$k_c$", y = "Counts") +
  theme_AP()

c <- ggplot(data.table(v.final), aes(v.final)) +
  geom_histogram() +
  labs(x = "$k_c$", y = "Counts") +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  theme_AP()

plot_grid(a, b, c, ncol = 3, labels = "auto")


## ----climate_uvalde, cache=TRUE--------------------------------------------------

# READ IN CLIMATIC DATA FOR UVALDE FOR JANUARY 2007-------------------------------

da <- fread("uvalde_climate_data_january.csv")

# Turn columns into numeric
col_names <- colnames(da)
da[, (col_names):= lapply(.SD, as.numeric), .SDcols = (col_names)]

# Convert Fahrenheit to Celsius
temp_cols <- colnames(da)[da[ , grepl( "Temp", colnames(da))]]
da[, (temp_cols):= lapply(.SD, function(x) (x - 32) / 1.8), .SDcols = (temp_cols)]

# Convert miles per hour to meters per second
wind_cols <- colnames(da)[da[ , grepl( "Wind", colnames(da))]]
da[, (wind_cols):= lapply(.SD, function(x) x / 2.237), .SDcols = (wind_cols)]

# Select days of interest
uvalde_days <- da[Day >= 6 & Day <= 7]

# Mean temperature
T_air <- mean(uvalde_days$`Temp.Avg`)
T_air

# Computation of Vapour pressure (Delta)
# T_air <- 27
e_sat <- 0.6108 * exp((17.27 * T_air) / (T_air + 237.3))
Delta <- (4098 * e_sat) / ((T_air + 237.3) ^ 2)
Delta

# Computation of vapour deficit (v)
rel_hum <- mean(uvalde_days$`Hum.Avg`) / 100
v <- e_sat - (e_sat * rel_hum)

# Computation of Psychrometric constant (gamma)
z <- 9.1
lambda <- 2.501 - 0.002361 * T_air
P <- 101.3 * ((293 - 0.0065 * z) / 293) ^ 5.256
Gamma <- 0.0016286 * P / lambda
Gamma


## ----settings_full, cache=TRUE---------------------------------------------------

# DEFINE SETTINGS ------------------------------------------------------------------

N <- 2^15
R <- 10^3
type <- "R"
order <- "second"
params <- c(
  "I_a", "Delta", "gamma", "A", "T_a", "w", "v",
  "k_c", "P", "E_a", "E_c", "M_f"
)

# Vector with the name of the parameters modified for better plotting
params.plot <- c(
  "$I_a$", "$\\Delta$", "$\\gamma$", "$A$", "$T_a$", "$w$",
  "$v$", "$k_c$", "$P$", "$E_a$", "$E_c$", "$M_f$"
)

# I_a (ha)
# P (mm)
# ET_0 (mm)


## ----sample_matrix_full, cache=TRUE, dependson=c("kc_wheat", "settings_full", "modify_rasters", "climate_uvalde")----

# DEFINE SAMPLE MATRIX AND TRANSFORM TO APPROPRIATE DISTRIBUTIONS -----------------

# Define sampling matrix
mat <- sobol_matrices(N = N, params = params, order = order, type = type)

# Transform to appropriate probability distributions
mat[, "I_a"] <- qunif(mat[, "I_a"], rasters.uvalde[, min(area)], rasters.uvalde[, max(area)])
mat[, "Delta"] <- qunif(mat[, "Delta"], Delta + Delta * -0.005, Delta + Delta * 0.005)
mat[, "gamma"] <- qunif(mat[, "gamma"], Gamma + Gamma * -0.001, Gamma + Gamma * 0.001)
mat[, "A"] <- qunif(mat[, "A"], 350 + 350 * -0.15, 350 + 350 * 0.15)
mat[, "T_a"] <- qunif(mat[, "T_a"], T_air + T_air * -0.01, T_air + T_air * 0.01)
mat[, "w"] <- qunif(mat[, "w"], 2.81 + 2.81 * -0.05, 2.81 + 2.81 * 0.05)
mat[, "v"] <- qunif(mat[, "v"], v + v * -0.04, v + v * 0.04)
mat[, "P"] <- qunif(mat[, "P"], 0, 0.1)
mat[, "M_f"] <- qunif(mat[, "M_f"], 0.5, 0.97)
mat[, "k_c"] <- qunif(mat[, "k_c"], min(kc_wheat.dt$y), max(kc_wheat.dt$y))
mat[, "E_a"] <- qunif(mat[, "E_a"], min(bos.dt[Type == "Sprinkler", `$E_a$`], na.rm = TRUE),
                      max(bos.dt[Type == "Sprinkler", `$E_a$`], na.rm = TRUE))
mat[, "E_c"] <- qunif(mat[, "E_c"], min(bos.dt[Type == "Sprinkler", `$E_c$`], na.rm = TRUE),
                      max(bos.dt[Type == "Sprinkler", `$E_c$`], na.rm = TRUE))


## ----plot_sample_matrix, cache=TRUE, dependson="sample_matrix_full", fig.height=6, fig.width=5.5----

# PLOT DISTRIBUTIONS ---------------------------------------------------------

data.table(mat[1:N, ]) %>%
  setnames(., params, params.plot) %>%
  melt(., measure.vars = params.plot) %>%
  ggplot(., aes(value)) +
  geom_histogram() +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  labs(x = "Value", y = "Counts") +
  theme_AP() +
  facet_wrap(~variable, scales = "free_x")


## ----define_model, cache=TRUE----------------------------------------------------

# DEFINE THE MODEL -----------------------------------------------------------

full_model <- function(I_a, Delta, A, gamma, T_a, w, v, k_c, P, E_a, E_c, M_f) {
  out <- ((k_c * ((0.408 * Delta * A + gamma * (900 / (T_a + 273)) * w * v) /
                    Delta + gamma * (1 + 0.34 * w))) - P) / (E_a * E_c * M_f)
  # Divide mm d-1 by 1000 to get m d-1
  y <- (out / 10^3) * I_a # Output is m3 ha
  return(y)
}


## ----run_model, cache=TRUE, dependson=c("define_model", "sample_matrix_full")----

# RUN THE MODEL -------------------------------------------------------------

y <- full_model(
  I_a = mat[, "I_a"],
  Delta = mat[, "Delta"],
  A = mat[, "A"],
  gamma = mat[, "gamma"],
  T_a = mat[, "T_a"],
  w = mat[, "w"],
  v = mat[, "v"],
  k_c = mat[, "k_c"],
  P = mat[, "P"],
  E_a = mat[, "E_a"],
  E_c = mat[, "E_c"],
  M_f = mat[, "M_f"]
)


## ----uncertainties_global, cache=TRUE, dependson="run_model"---------------------

# ASSESS UNCERTAINTIES -------------------------------------------------------

unc <- plot_uncertainty(Y = y, N = N)


## ----sensitivities_full, cache=TRUE, dependson="run_model"-----------------------

# ASSESS SENSITIVITIES ------------------------------------------------------------

# Compute sobol' indices
ind <- sobol_indices(
  Y = y, N = N, params = params.plot,
  order = order, boot = TRUE, R = R,
  parallel = "multicore"
)

# Plot
ind

# Everything is explained by first and second-order effects
ind$results[sensitivity %in% c("Si", "Sij"), sum(original)]

# Plot sobol' indices
sobol.plot <- plot(ind) +
  theme(legend.position = c(0.83, 0.5))


## ----second_order, cache=TRUE, dependson="sensitivities_full", fig.width=4, fig.height=2.4----

# PLOT SECOND-ORDER INDICES -------------------------------------------------------

second.order <- plot(ind, "second")


## ----plot_sensitivities_scatter, cache=TRUE, dependson=c("sample_matrix_full", "settings_full", "run_model"), fig.height=6, fig.width=5.5----

# PLOT SCATTERPLOTS ---------------------------------------------------------------


## ----oat_matrix, cache=TRUE, dependson=c("sample_matrix_full", "run_model")------

# CONSTRUCT SAMPLE MATRIX ---------------------------------------------------------

A <- mat[1:N, ]
B <- matrix(rep(Rfast::colmeans(A), each = N), nrow = N)

X <- B
for (j in 1:ncol(A)) {
  AB <- B
  AB[, j] <- A[, j]
  X <- rbind(X, AB)
}

mat.oat <- X[(N + 1):nrow(X), ]
colnames(mat.oat) <- params


## ----model_oat, cache=TRUE, dependson=c("oat_matrix", "run_model")---------------

# RUN THE MODEL ------------------------------------------------------------------

y.oat <- full_model(
  I_a = mat.oat[, "I_a"],
  Delta = mat.oat[, "Delta"],
  A = mat.oat[, "A"],
  gamma = mat.oat[, "gamma"],
  T_a = mat.oat[, "T_a"],
  w = mat.oat[, "w"],
  v = mat.oat[, "v"],
  k_c = mat.oat[, "k_c"],
  P = mat.oat[, "P"],
  E_a = mat.oat[, "E_a"],
  E_c = mat.oat[, "E_c"],
  M_f = mat.oat[, "M_f"]
)


## ----single_point, cache=TRUE, dependson=c("sample_matrix_full", "oat_matrix")----

# COMPUTE A SINGLE-POINT ESTIMATE USING MEAN VALUES-------------------------------

vec_means <- colMeans(A)

y.point <- full_model(
  I_a = vec_means[["I_a"]],
  Delta = vec_means[["Delta"]],
  A = vec_means[["A"]],
  gamma = vec_means[["gamma"]],
  T_a = vec_means[["T_a"]],
  w = vec_means[["w"]],
  v = vec_means[["v"]],
  k_c = vec_means[["k_c"]],
  P = vec_means[["P"]],
  E_a = vec_means[["E_a"]],
  E_c = vec_means[["E_c"]],
  M_f = vec_means[["M_f"]]
)

y.point


## ----unc_oat, cache=TRUE, dependson=c("model_oat", "settings_full", "sample_matrix_full"), fig.height=3, fig.width=4.5----

# ASSESS UNCERTAINTIES --------------------------------------------------------

unc.oat <- plot_uncertainty(Y = y.oat, N = N)

full.unc <- data.table(cbind(y[1:N], y.oat))
colnames(full.unc) <- c("Global", "OAT")

a <- melt(full.unc, measure.vars = colnames(full.unc)) %>%
  ggplot(., aes(value, fill = variable)) +
  geom_histogram(position = "identity", alpha = 0.3, color = "black") +
  labs(x = "Irrigation water withdrawal (m$^3$ ha d$^{-1}$)", y = "Counts") +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_vline(xintercept = y.point, lty = 2, color = "red", size = 2) +
  scale_fill_manual(values = wes_palette(2, name = "Chevalier1"), name = "Uncertainty analysis") +
  theme_AP() +
  theme(legend.position = c(0.2, 0.5))

a

d <- melt(full.unc, measure.vars = colnames(full.unc)) %>%
  .[variable == "Global"] %>%
  ggplot(., aes(value)) +
  geom_histogram(position = "identity", alpha = 0.3, color = "black", fill = "white") +
  labs(x = "Irrigation water withdrawal (m$^3$ ha d$^{-1}$)", y = "Counts") +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_vline(xintercept = y.point, lty = 2, color = "red", size = 2) +
  scale_fill_manual(values = wes_palette(2, name = "Chevalier1"), name = "Uncertainty analysis") +
  theme_AP() +
  theme(legend.position = c(0.2, 0.5))

d

all.projection <- d # For later


## ----stats_unc, cache=TRUE, dependson="unc_oat"----------------------------------

# SOME STATISTICS --------------------------------------------------------------

stat.full.unc <- melt(full.unc, measure.vars = c("Global", "OAT"))
stat.full.unc[, .(min = min(value), max = max(value)), variable]

# Quantiles
probs.quantile <- c(0, 0.025, 0.1, 0.5, 0.9, 0.975, 1)
stat.full.unc[, .(value = quantile(value, probs = probs.quantile)), variable] %>%
  .[, quantile := rep(probs.quantile, 2)] %>%
  dcast(., variable ~ quantile, value.var = "value") %>%
  print()


## ----merge_full_oat, cache=TRUE, dependson=c("uncertainties_global", "unc_oat", "sensitivities_full", "second_order"), fig.height = 6, fig.width = 4.7----

# MERGE UNCERTAINTY AND SOBOL' INDICES ---------------------------------------

plot_grid(a, sobol.plot, second.order, ncol = 1, labels = "auto")


## ----all_figures_comment, cache=TRUE, fig.height=4.5, fig.width=4.5--------------

# PLOT --------------------------------------------------------------------------

legend <- get_legend(plots.unc[[2]] + theme(legend.position = "top"))
bottom <- plot_grid(etc + theme(legend.position = "none"), ep, labels = "auto")
top <- plot_grid(legend, bottom, rel_heights = c(0.13, 0.87), ncol = 1)
plot_grid(top, all.projection, labels = c("", "c"), ncol = 1)


