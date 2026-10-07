# Custom functions for the vocal character displacement analyses.
# Sourced by every report in scripts/ (see the "load packages" chunk).

print <- function(x, row.names = FALSE) {
  kb <- kable(x, row.names = row.names, digits = 4, "html")
  kb <- kable_styling(kb,
                      bootstrap_options = c("striped", "hover", "condensed", "responsive"))
  scroll_box(kb, width = "100%")
}

get_stable_loadings <- function(pca, pcs = 4, B = 1000, cum_threshold = 0.5,
    freq_threshold = 0.75, seed = 123) {

    set.seed(seed)

    # Reconstruct centered data
    X <- pca$x %*% t(pca$rotation)

    p <- ncol(pca$rotation)
    orig_rot <- pca$rotation[, 1:pcs]

    boot_loadings <- array(NA, dim = c(p, pcs, B))

    # ----------------------------- Bootstrap PCA
    # -----------------------------

    for (b in 1:B) {

        idx <- sample(1:nrow(X), replace = TRUE)
        Xb <- X[idx, ]

        pca_b <- prcomp(Xb, scale. = TRUE)
        rot_b <- pca_b$rotation[, 1:pcs]

        # Align signs
        for (k in 1:pcs) {
            if (cor(rot_b[, k], orig_rot[, k]) < 0) {
                rot_b[, k] <- -rot_b[, k]
            }
        }

        boot_loadings[, , b] <- rot_b
    }

    # ----------------------------- Summary statistics
    # -----------------------------

    abs_boot <- abs(boot_loadings)

    mean_loading <- apply(abs_boot, c(1, 2), mean)
    ci_lower <- apply(abs_boot, c(1, 2), quantile, 0.025)
    ci_upper <- apply(abs_boot, c(1, 2), quantile, 0.975)

    # ----------------------------- Stability frequency
    # -----------------------------

    top_freq <- matrix(0, nrow = p, ncol = pcs)

    for (b in 1:B) {
        for (k in 1:pcs) {

            sq <- boot_loadings[, k, b]^2
            ord <- order(sq, decreasing = TRUE)
            cumprop <- cumsum(sq[ord])/sum(sq)

            selected <- ord[cumprop <= cum_threshold]
            selected <- c(selected, ord[min(which(cumprop >= cum_threshold))])

            top_freq[selected, k] <- top_freq[selected, k] + 1
        }
    }

    top_freq <- top_freq/B

    stable <- (top_freq >= freq_threshold) & (ci_lower > 0 | ci_upper <
        0)

    # ----------------------------- Return tidy dataframe
    # -----------------------------

    data.frame(variable = rep(rownames(pca$rotation), pcs), ind = rep(colnames(pca$rotation)[1:pcs],
        each = p), mean_loading = as.vector(mean_loading), ci_lower = as.vector(ci_lower),
        ci_upper = as.vector(ci_upper), freq = as.vector(top_freq),
        stable = as.vector(stable))
}

# number of PCs to retain using the broken-stick criterion (Frontier 1976; Jackson 1993):
# a PC is retained while its observed proportion of variance exceeds the proportion
# expected under a random division ("broken stick") of total variance among the same
# number of components
n_pcs_broken_stick <- function(pca) {

    # observed proportion of variance explained by each PC
    obs_var <- summary(pca)$importance[2, ]

    p <- length(obs_var)

    # expected proportion of variance under the broken-stick null model
    bstick_expected <- sapply(seq_len(p), function(k) sum(1 / (k:p)) / p)

    # first PC (if any) at which the observed variance no longer exceeds
    # the broken-stick expectation
    below <- which(obs_var <= bstick_expected)

    n_keep <- if (length(below) == 0) p else below[1] - 1

    # always retain at least one PC
    max(n_keep, 1)
}

# higliht significant rows

highlight <- function(x, estimate_col = "Estimate", lower_col = "Q2.5",
    upper_col = "Q97.5", strong_fill = "#E8602DFF", moderate_fill = "#FAC127FF",
    weak_fill = "#FCFFA4FF", alpha = 0.5, digits = 3) {
    ## -------------------------------------------------- Row
    ## groups --------------------------------------------------

    strong_rows <- which(x$pd > 0.95)

    moderate_rows <- which(x$pd > 0.9 & x$pd <= 0.95)

    weak_rows <- which(x$pd > 0.8 & x$pd <= 0.9)

    ## -------------------------------------------------- Build
    ## kable --------------------------------------------------

    x_kbl <- kableExtra::kbl(x, row.names = TRUE, escape = FALSE,
        format = "html", digits = digits)

    ## -------------------------------------------------- Apply
    ## row highlighting
    ## --------------------------------------------------

    if (length(strong_rows) > 0) {

        x_kbl <- kableExtra::row_spec(x_kbl, row = strong_rows, background = grDevices::adjustcolor(strong_fill,
            alpha.f = alpha))
    }

    if (length(moderate_rows) > 0) {

        x_kbl <- kableExtra::row_spec(x_kbl, row = moderate_rows,
            background = grDevices::adjustcolor(moderate_fill, alpha.f = alpha))
    }

    if (length(weak_rows) > 0) {

        x_kbl <- kableExtra::row_spec(x_kbl, row = weak_rows, background = grDevices::adjustcolor(weak_fill,
            alpha.f = alpha))
    }

    ## --------------------------------------------------
    ## Styling
    ## --------------------------------------------------

    x_kbl <- kableExtra::kable_styling(x_kbl, bootstrap_options = c("striped",
        "hover", "condensed", "responsive"), full_width = FALSE, font_size = 12)

    return(x_kbl)
}

plot_brms_heatmap <- function(
    model_files,
    remove_intercepts = TRUE
) {

  # models may still be fitting in the background, so some (or all) of
  # the expected files can be missing, and a file that is still being
  # written by brms can be present but not yet readable - skip either
  # case instead of erroring
  if(length(model_files) == 0) {
    message("plot_brms_heatmap: no model files found yet - skipping.")
    return(invisible(NULL))
  }

  effects_df <- data.frame()

  for(i in seq_along(model_files)) {

    fit <- tryCatch(
      readRDS(model_files[i]),
      error = function(e) NULL
    )

    if(is.null(fit)) {
      message(
        "plot_brms_heatmap: could not read '", model_files[i],
        "' (likely still being written) - skipping."
      )
      next
    }

    fe <- as.data.frame(fixef(fit))
    fe$predictor <- rownames(fe)

    if(remove_intercepts)
      fe <- fe[fe$predictor != "Intercept", ]

response <- deparse(fit$formula$formula[[2]])

    fe$response <- response

    effects_df <- rbind(
      effects_df,
      fe
    )
  }

  if(nrow(effects_df) == 0) {
    message("plot_brms_heatmap: no readable model results yet - skipping.")
    return(invisible(NULL))
  }

  rownames(effects_df) <- NULL

  # significance

  effects_df$sig <- with(
    effects_df,
    Q2.5 * Q97.5 > 0
  )

  # clean names

  effects_df$predictor <- gsub("\\)$", "", effects_df$predictor)

  effects_df$predictor <- gsub("^mo", "", effects_df$predictor)
  effects_df$predictor <- gsub("^mi", "", effects_df$predictor)
  effects_df$predictor <- gsub("_sc$", "", effects_df$predictor)

  effects_df$response <- gsub("^mi", "", effects_df$response)
  effects_df$response <- gsub("_sc$", "", effects_df$response)

  # average duplicated cells if present

  effects_df <- aggregate(
    cbind(
      Estimate,
      sig
    ) ~ predictor + response,
    data = effects_df,
    FUN = mean
  )

  effects_df$sig <- effects_df$sig > 0.5
  
  # complete combinations

  all_combos <- expand.grid(
    predictor = unique(effects_df$predictor),
    response = unique(effects_df$response)
  )

  plot_df <- merge(
    all_combos,
    effects_df,
    by = c("predictor", "response"),
    all.x = TRUE
  )

 plot_df$response <- gsub("acoustic_distance$", "acoustic\ndistance", plot_df$response)

 plot_df$response <- gsub("_distance$", "", plot_df$response)

plot_df$predictor <- gsub("sympatry1", "Sympatry", plot_df$predictor)
    
  lim <- max(abs(plot_df$Estimate), na.rm = TRUE)

  
  ggplot(
    plot_df,
    aes(
      predictor,
      response
    )
  ) +

    geom_tile(
      fill = "grey90",
      colour = "white"
    ) +

    geom_tile(
      data = plot_df[!is.na(plot_df$Estimate), ],
      aes(fill = Estimate),
      colour = "white"
    ) +

    geom_text(
      data = plot_df[!is.na(plot_df$Estimate), ],
      aes(
        label = sprintf("%.2f", Estimate),
        colour = sig
      ),
      fontface = "bold",
      size = 3
    ) +
      
    # scale_fill_gradient2(
    #   low = rep("#403B78", 2),
    #   mid = "white",
    #   high = rep("#DEF5E5", rep = 2),
    #   midpoint = 0,
    #   name = "Estimate"
    # ) +
    # 
      scale_fill_gradient2(
  low = "#403B78",
  mid = "white",
  high = "#A0DFB9CC",
  midpoint = 0,
  limits = c(-lim, lim),
  oob = scales::squish,
  name = "Estimate"
) +

    scale_color_manual(
      values = c(
        "TRUE" = "black",
        "FALSE" = "grey70"
      ),
      guide = "none"
    ) +

    labs(
      x = "Predictor",
      y = "Response"
    ) +

    theme_classic() +

    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
}

# table of fixed-effect estimates (posterior mean, SE and 95% uncertainty
# interval) for a set of saved brms fits, one row per model x predictor;
# skips missing/unreadable files so it can be run while models are still
# fitting
fixef_table <- function(model_files, remove_intercepts = TRUE, digits = 3) {

  rows <- lapply(model_files, function(f) {
    fit <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    fe <- as.data.frame(fixef(fit))
    fe$predictor <- rownames(fe)
    if (remove_intercepts) fe <- fe[fe$predictor != "Intercept", ]
    fe$response <- deparse(fit$formula$formula[[2]])
    rm(fit)
    invisible(gc(verbose = FALSE))
    fe
  })

  out <- do.call(rbind, rows)

  if (is.null(out) || nrow(out) == 0) {
    message("fixef_table: no readable model results yet - skipping.")
    return(invisible(NULL))
  }

  rownames(out) <- NULL

  out$response  <- gsub("_distance(_sc)?$", "", out$response)
  out$predictor <- gsub("^sympatry1$", "Sympatry", out$predictor)
  out$ui_excludes_zero <- ifelse(out$Q2.5 > 0 | out$Q97.5 < 0, "yes", "no")

  out <- out[, c("response", "predictor", "Estimate", "Est.Error", "Q2.5", "Q97.5", "ui_excludes_zero")]
  out[, 3:6] <- round(out[, 3:6], digits)
  names(out) <- c("Response", "Predictor", "Estimate", "SE", "l-95% UI", "u-95% UI", "UI excludes 0")

  out[order(out$Response, out$Predictor), ]
}

# ---- divergent transitions diagnostics ------------------------------------
# where do the divergent transitions of a saved brms fit occur, and do the
# fixed-effect estimates depend on them? Uses the existing fit only (no
# refitting). Returns a list with:
#   summary : number and % of divergent post-warmup draws, divergences per
#             chain and number of draws that hit the maximum treedepth
#   fixef   : 2.5%, 50% and 97.5% posterior quantiles of each fixed effect,
#             separately for divergent and non-divergent draws
#   sd      : median of each random-effect SD in non-divergent vs divergent
#             draws, plus the mean percentile of the divergent draws within
#             the whole posterior of that SD (~0.5 = divergences spread over
#             the posterior; close to 0 = divergences piled up at small SD
#             values, i.e. a funnel involving that SD)
# If plot = TRUE, also draws each random-effect SD (log scale) against the
# focal coefficient with divergent draws in red.
divergence_check <- function(fit, focal = "b_sympatry1", max_treedepth = 15,
                             digits = 3, plot = TRUE) {

  pars <- grep("^b_|^sd_", variables(fit), value = TRUE)

  draws <- as.data.frame(as_draws_df(fit, variable = pars))

  np <- nuts_params(fit)

  div_np <- np[np$Parameter == "divergent__", c("Chain", "Iteration", "Value")]
  names(div_np) <- c(".chain", ".iteration", "divergent")

  td_np <- np[np$Parameter == "treedepth__", c("Chain", "Iteration", "Value")]
  names(td_np) <- c(".chain", ".iteration", "treedepth")

  draws <- merge(draws, div_np, by = c(".chain", ".iteration"))
  draws <- merge(draws, td_np, by = c(".chain", ".iteration"))

  div <- draws$divergent == 1

  # overall summary
  per_chain <- tapply(draws$divergent, draws$.chain, sum)

  summ <- data.frame(
    `divergent draws` = sum(div),
    `total draws` = nrow(draws),
    `% divergent` = round(100 * mean(div), 2),
    `divergent per chain` = paste(per_chain, collapse = " / "),
    `max treedepth hits` = sum(draws$treedepth >= max_treedepth),
    check.names = FALSE
  )

  # fixed effects in divergent vs non-divergent draws
  q <- function(x) {
    if (length(x) == 0) return(c(`2.5%` = NA, `50%` = NA, `97.5%` = NA))
    quantile(x, c(0.025, 0.5, 0.975))
  }

  b_pars <- setdiff(grep("^b_", pars, value = TRUE), "b_Intercept")

  fe <- do.call(rbind, lapply(b_pars, function(p) {
    data.frame(
      parameter = p,
      draws = c("non-divergent", "divergent"),
      n = c(sum(!div), sum(div)),
      rbind(q(draws[[p]][!div]), q(draws[[p]][div])),
      check.names = FALSE
    )
  }))
  fe[, 4:6] <- round(fe[, 4:6], digits)

  # where do divergences fall along each random-effect SD?
  sd_pars <- grep("^sd_", pars, value = TRUE)

  sds <- do.call(rbind, lapply(sd_pars, function(p) {
    x <- draws[[p]]
    data.frame(
      parameter = p,
      `median (non-divergent)` = median(x[!div]),
      `median (divergent)` = if (any(div)) median(x[div]) else NA,
      `mean percentile of divergent draws` = if (any(div)) mean(ecdf(x)(x[div])) else NA,
      check.names = FALSE
    )
  }))
  sds[, 2:4] <- round(sds[, 2:4], digits)

  # scatter plots: each SD (log) vs the focal coefficient
  if (plot && any(div) && focal %in% pars) {

    arr <- as.array(fit, variable = c(sd_pars, focal))

    for (p in sd_pars) {
      # base::print - the document redefines print() as a kable wrapper
      base::print(
        bayesplot::mcmc_scatter(
          arr,
          pars = c(p, focal),
          np = np,
          transformations = setNames(list("log"), p),
          alpha = 0.15,
          size = 0.6,
          np_style = bayesplot::scatter_style_np(div_color = "red", div_size = 1.2, div_alpha = 0.8)
        ) +
          labs(title = paste("Divergent draws (red):", p), x = paste0("log(", p, ")"), y = focal) +
          theme_classic(base_size = 12)
      )
    }
  }

  invisible(list(summary = summ, fixef = fe, sd = sds))
}

# divergence summary for a set of saved fits (one row per model), with the
# focal coefficient's median in non-divergent vs divergent draws; skips
# missing/unreadable files
divergence_table <- function(model_files, focal = "b_sympatry1", max_treedepth = 15, digits = 3) {

  rows <- lapply(model_files, function(f) {
    fit <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    dc <- divergence_check(fit, focal = focal, max_treedepth = max_treedepth,
                           digits = digits, plot = FALSE)

    foc <- dc$fixef[dc$fixef$parameter == focal, ]
    sd_sp <- dc$sd[grepl("species_pair", dc$sd$parameter), ]

    out <- data.frame(
      response = gsub("_distance(_sc)?$", "", deparse(fit$formula$formula[[2]])),
      dc$summary,
      `sympatry median (non-divergent)` = foc$`50%`[foc$draws == "non-divergent"],
      `sympatry median (divergent)` = foc$`50%`[foc$draws == "divergent"],
      `species-pair SD: mean percentile of divergent draws` =
        if (nrow(sd_sp)) sd_sp$`mean percentile of divergent draws`[1] else NA,
      check.names = FALSE
    )

    rm(fit)
    invisible(gc(verbose = FALSE))
    out
  })

  out <- do.call(rbind, rows)

  if (is.null(out) || nrow(out) == 0) {
    message("divergence_table: no readable model results yet - skipping.")
    return(invisible(NULL))
  }

  out[order(out$response), ]
}

# ---- contrasts by species pair -------------------------------------------
# sympatry effect for each species pair from a model with species-pair-
# specific sympatry slopes, i.e. (1 + sympatry | species_pair): per-pair
# effect = fixed sympatry effect + the pair's deviation. Returns the per-pair
# effects (plus the pooled fixed effect), the pairwise differences between
# species pairs, and a plot.
species_pair_sympatry_effects <- function(fit, digits = 3) {

  # summary per species pair (fixed + random slope)
  cf <- coef(fit, summary = TRUE)$species_pair[, , "sympatry1", drop = TRUE]
  cf <- as.data.frame(cf)
  cf$species_pair <- rownames(cf)
  rownames(cf) <- NULL

  pooled <- fixef(fit)["sympatry1", ]
  cf <- rbind(
    cf,
    data.frame(
      Estimate = pooled[["Estimate"]], Est.Error = pooled[["Est.Error"]],
      Q2.5 = pooled[["Q2.5"]], Q97.5 = pooled[["Q97.5"]],
      species_pair = "Pooled (fixed effect)"
    )
  )
  cf$ui_excludes_zero <- ifelse(cf$Q2.5 > 0 | cf$Q97.5 < 0, "yes", "no")
  effects <- cf[, c("species_pair", "Estimate", "Est.Error", "Q2.5", "Q97.5", "ui_excludes_zero")]
  effects[, 2:5] <- round(effects[, 2:5], digits)
  names(effects) <- c("Species pair", "Sympatry effect", "SE", "l-95% UI", "u-95% UI", "UI excludes 0")

  # pairwise differences between species pairs, from the posterior draws
  draws <- coef(fit, summary = FALSE)$species_pair[, , "sympatry1"]
  lev <- colnames(draws)
  diffs <- list()
  if (length(lev) > 1) {
    cmb <- combn(lev, 2)
    for (k in seq_len(ncol(cmb))) {
      d <- draws[, cmb[1, k]] - draws[, cmb[2, k]]
      diffs[[k]] <- data.frame(
        contrast = paste(cmb[1, k], "-", cmb[2, k]),
        Estimate = mean(d),
        Q2.5 = unname(quantile(d, 0.025)),
        Q97.5 = unname(quantile(d, 0.975))
      )
    }
  }
  diffs <- do.call(rbind, diffs)
  if (!is.null(diffs)) {
    diffs$ui_excludes_zero <- ifelse(diffs$Q2.5 > 0 | diffs$Q97.5 < 0, "yes", "no")
    diffs[, 2:4] <- round(diffs[, 2:4], digits)
    names(diffs) <- c("Contrast", "Difference", "l-95% UI", "u-95% UI", "UI excludes 0")
  }

  plot_df <- cf
  plot_df$species_pair <- factor(plot_df$species_pair, levels = rev(cf$species_pair))
  plot_df$type <- ifelse(plot_df$species_pair == "Pooled (fixed effect)", "Pooled", "Species pair")

  p <- ggplot(plot_df, aes(x = Estimate, y = species_pair, color = type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5), width = 0, orientation = "y", linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Species pair" = "#403B78", "Pooled" = "#5FA98A"), guide = "none") +
    labs(x = "Sympatry effect (SD of acoustic distance) with 95% UI", y = NULL) +
    theme_classic()

  list(effects = effects, differences = diffs, plot = p)
}

# ---- contrasts among population pairs -------------------------------------
# expected (standardized) acoustic distance for every observed population-
# pair combination within each species pair, from the fitted model
# including the species-pair and population (multi-membership) effects but
# not the individual effects; then sympatric vs allopatric population-pair
# contrasts within species pairs.
population_pair_contrasts <- function(fit, dat, ndraws = 2000, digits = 3,
                                      show_posterior = TRUE) {

  # put the two populations of each comparison in a fixed (alphabetical)
  # order, so that A-B and B-A comparisons are pooled into one population
  # pair (the multiple-membership term is symmetric, so order does not
  # affect the prediction)
  dat <- as.data.frame(dat)
  p1 <- as.character(dat$pop1)
  p2 <- as.character(dat$pop2)
  swap <- p1 > p2
  dat$pop_a <- ifelse(swap, p2, p1)
  dat$pop_b <- ifelse(swap, p1, p2)
  dat$n_comparisons <- 1

  nd <- aggregate(
    n_comparisons ~ species_pair + pop_a + pop_b + sympatry,
    data = dat, FUN = sum
  )
  names(nd)[names(nd) == "pop_a"] <- "pop1"
  names(nd)[names(nd) == "pop_b"] <- "pop2"

  # individual columns are required by brms but excluded from the prediction
  nd$individual1 <- dat$individual1[1]
  nd$individual2 <- dat$individual2[1]

  # keep the species-pair and population terms, drop the individual terms;
  # if the model has species-pair-specific sympatry slopes, keep them too
  has_sp_slopes <- "sympatry1" %in% dimnames(ranef(fit)$species_pair)[[3]]
  re_form <- if (has_sp_slopes) {
    ~ (1 + sympatry | species_pair) + (1 | mm(pop1, pop2))
  } else {
    ~ (1 | species_pair) + (1 | mm(pop1, pop2))
  }

  ep <- posterior_epred(
    fit,
    newdata = nd,
    re_formula = re_form,
    ndraws = ndraws
  )

  nd$Estimate <- colMeans(ep)
  nd$Q2.5 <- apply(ep, 2, quantile, 0.025)
  nd$Q97.5 <- apply(ep, 2, quantile, 0.975)

  # pop1/pop2 are "species.population" (interaction(), needed to keep
  # population codes that repeat across species unique for the mm() term);
  # for display, drop the species prefix - the facet strip (species_pair)
  # already identifies which two species are being compared
  strip_species <- function(x) sub("^[^.]*\\.", "", as.character(x))
  nd$population_pair <- paste(strip_species(nd$pop1), strip_species(nd$pop2), sep = " vs ")
  nd$sympatry_label <- ifelse(nd$sympatry == "1", "Sympatric", "Allopatric")

  # sympatric vs allopatric population pairs within each species pair
  contrasts <- list()
  for (sp in unique(as.character(nd$species_pair))) {
    idx_s <- which(nd$species_pair == sp & nd$sympatry == "1")
    idx_a <- which(nd$species_pair == sp & nd$sympatry == "0")
    for (i in idx_s) for (j in idx_a) {
      d <- ep[, i] - ep[, j]
      contrasts[[length(contrasts) + 1]] <- data.frame(
        species_pair = sp,
        sympatric_pair = nd$population_pair[i],
        allopatric_pair = nd$population_pair[j],
        Estimate = mean(d),
        Q2.5 = unname(quantile(d, 0.025)),
        Q97.5 = unname(quantile(d, 0.975))
      )
    }
  }
  contrasts <- do.call(rbind, contrasts)
  if (!is.null(contrasts)) {
    contrasts$ui_excludes_zero <- ifelse(contrasts$Q2.5 > 0 | contrasts$Q97.5 < 0, "yes", "no")
    contrasts[, 4:6] <- round(contrasts[, 4:6], digits)
    names(contrasts) <- c("Species pair", "Sympatric population pair", "Allopatric population pair",
                          "Difference", "l-95% UI", "u-95% UI", "UI excludes 0")
  }

  expected <- nd[, c("species_pair", "population_pair", "sympatry_label", "n_comparisons",
                     "Estimate", "Q2.5", "Q97.5")]
  expected[, 5:7] <- round(expected[, 5:7], digits)
  names(expected) <- c("Species pair", "Population pair", "Sympatry", "N comparisons",
                       "Expected distance", "l-95% UI", "u-95% UI")

  # one row of nd per (species_pair, population pair) combination, so the
  # count per species_pair (in its factor-level order, matching the order
  # facet_wrap2() draws panels in) is exactly the number of population-pair
  # rows that panel needs
  panel_counts <- as.numeric(table(nd$species_pair))

  # ggh4x::facet_wrap2() keeps the facet_wrap()-style banner strip across
  # the top of each panel (facet_grid()'s "space = free" doesn't), while
  # ggh4x::force_panelsizes() still gives each panel a height proportional
  # to how many population-pair rows it holds, instead of every panel
  # getting the same height regardless of row count
  # population pairs ordered by their mean expected distance (same order
  # reorder(population_pair, Estimate) gave), shared by the point/interval
  # layer and the posterior slabs so both land on the same rows
  pair_levels <- names(sort(tapply(nd$Estimate, nd$population_pair, mean)))
  nd$population_pair_f <- factor(nd$population_pair, levels = pair_levels)

  sym_cols <- c("Allopatric" = "#403B78", "Sympatric" = "#5FA98A")

  p <- ggplot(nd, aes(x = Estimate, y = population_pair_f, color = sympatry_label))

  # posterior distribution of the expected distance for each population pair
  # (one half-density per row, drawn behind the mean and 95% UI); each slab
  # is scaled to its own maximum so that narrow and wide posteriors are both
  # visible - slab height is therefore not comparable across rows, only shape
  # and spread are
  if (show_posterior) {
    draws_long <- data.frame(
      species_pair = rep(nd$species_pair, each = nrow(ep)),
      population_pair_f = factor(rep(nd$population_pair, each = nrow(ep)), levels = pair_levels),
      sympatry_label = rep(nd$sympatry_label, each = nrow(ep)),
      value = as.vector(ep)  # column-major: all draws of row 1, then row 2, ...
    )

    p <- p +
      ggdist::stat_slab(
        data = draws_long,
        aes(x = value, y = population_pair_f, fill = sympatry_label),
        inherit.aes = FALSE,
        normalize = "groups",
        scale = 0.8,
        alpha = 0.35,
        color = NA
      ) +
      scale_fill_manual(values = sym_cols, guide = "none")
  }

  p <- p +
    geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5), width = 0, orientation = "y", linewidth = 1) +
    geom_point(size = 3) +
    ggh4x::facet_wrap2(~ species_pair, scales = "free_y", ncol = 1) +
    ggh4x::force_panelsizes(rows = panel_counts) +
    scale_color_manual(values = sym_cols, name = NULL) +
    labs(
      x = if (show_posterior) "Expected acoustic distance (SD): posterior, mean and 95% UI"
          else "Expected acoustic distance (SD) with 95% UI",
      y = NULL
    ) +
    theme_classic() +
    theme(legend.position = "top")

  list(expected = expected, contrasts = contrasts, plot = p)
}

# ---- population-pair contrasts across every per-feature model -------------
# runs population_pair_contrasts() once per saved per-feature fit (the same
# model_files list assembled for plot_brms_heatmap()/fixef_table() next to
# the feature-based heatmap), tags each resulting contrast with the feature
# it came from, and returns the combined contrasts table plus a named list
# of plots - one per feature, not one plot faceted across features - each
# reusing population_pair_contrasts()'s own plot (which already facets by
# species pair) so every feature gets its own, fully legible figure. Skips
# any file that fails to load or fails to fit (e.g. a model still running,
# or one for which population_pair_contrasts() errors) rather than
# stopping the whole loop.
population_pair_contrasts_by_feature <- function(model_files, dat, ndraws = 2000, digits = 3) {

  all_contrasts <- list()
  plots <- list()

  for (f in model_files) {

    if (!file.exists(f)) next

    feature_name <- tools::file_path_sans_ext(basename(f))
    feature_name <- gsub("_distance(_sc)?_sympatry_only.*", "", feature_name)
    feature_name <- gsub("_", " ", feature_name)

    fit <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(fit)) next

    res <- tryCatch(
      population_pair_contrasts(fit = fit, dat = dat, ndraws = ndraws, digits = digits),
      error = function(e) NULL
    )

    rm(fit)

    if (is.null(res) || is.null(res$contrasts)) next

    res$contrasts$Feature <- feature_name
    all_contrasts[[feature_name]] <- res$contrasts
    plots[[feature_name]] <- res$plot + labs(title = feature_name)
  }

  invisible(gc(verbose = FALSE))

  if (length(all_contrasts) == 0) return(NULL)

  contrasts <- do.call(rbind, all_contrasts)
  rownames(contrasts) <- NULL
  contrasts <- contrasts[, c("Feature", "Species pair", "Sympatric population pair",
                              "Allopatric population pair", "Difference", "l-95% UI",
                              "u-95% UI", "UI excludes 0")]

  list(contrasts = contrasts, plots = plots)
}

# correlation check among the song-level features that feed the PCA and
# are used as responses in the per-feature models: pairwise Pearson
# correlations among songs. Prints the pairs above the correlation
# threshold, draws a correlation matrix plot, and returns the matrix and
# pair table invisibly.
check_feature_collinearity <- function(
    dat,
    features,
    labels = NULL,
    r_threshold = 0.7
) {

  X <- dat[, features, drop = FALSE]
  X <- X[complete.cases(X), ]
  X[] <- lapply(X, as.numeric)

  if (is.null(labels)) labels <- features
  names(labels) <- features

  # pairwise correlations
  cor_mat <- cor(X, method = "pearson")

  cor_long <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  names(cor_long) <- c("feature_1", "feature_2", "r")
  cor_long$feature_1 <- factor(labels[cor_long$feature_1], levels = labels)
  cor_long$feature_2 <- factor(labels[cor_long$feature_2], levels = labels)

  # unique pairs above threshold
  pair_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
  pairs_df <- data.frame(
    feature_1 = labels[rownames(cor_mat)[pair_idx[, 1]]],
    feature_2 = labels[colnames(cor_mat)[pair_idx[, 2]]],
    r = round(cor_mat[pair_idx], 3),
    row.names = NULL
  )
  pairs_df <- pairs_df[order(-abs(pairs_df$r)), ]
  high_pairs <- pairs_df[abs(pairs_df$r) >= r_threshold, ]

  cat("\nPairs with |r| >=", r_threshold, ":\n")
  if (nrow(high_pairs) == 0) {
    cat("  none\n")
  } else {
    base::print(high_pairs, row.names = FALSE)
  }

  cat("\nStrongest pairwise correlations:\n")
  base::print(head(pairs_df, 5), row.names = FALSE)

  # correlation matrix plot (lower triangle), same palette as the model heatmaps
  cor_long$show <- as.integer(cor_long$feature_1) > as.integer(cor_long$feature_2)
  cor_long$flag <- abs(cor_long$r) >= r_threshold

  p <- ggplot(cor_long[cor_long$show, ], aes(x = feature_2, y = feature_1, fill = r)) +
    geom_tile(color = "white") +
    geom_text(
      aes(
        label = sprintf("%.2f", r),
        fontface = ifelse(flag, "bold", "plain"),
        color = flag
      ),
      size = 3.2
    ) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey40"), guide = "none") +
    scale_fill_gradient2(
      low = "#403B78",
      mid = "white",
      high = "#A0DFB9CC",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson r"
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # base::print - the document redefines print() as a kable wrapper
  base::print(p)

  invisible(list(cor = cor_mat, pairs = pairs_df, high_pairs = high_pairs))
}

# ===========================================================================
# within-species variation (Steps 1 and 3)
# ===========================================================================

# great-circle distance (km) between points given in decimal degrees
haversine_km <- function(lat1, lon1, lat2, lon2) {
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  6371 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

# every unique song x song pair (conspecific and heterospecific) with the
# PCA-based Euclidean acoustic distance (same definition as in the sympatry
# models), species / population / individual of each song, species-specific
# population IDs and geographic distance between the two recordings
all_song_pairs <- function(df, pcs) {

  d <- as.matrix(dist(df[, pcs], method = "euclidean"))
  idx <- which(upper.tri(d), arr.ind = TRUE)
  i <- idx[, 1]
  j <- idx[, 2]

  out <- data.frame(
    id1 = as.character(df$song[i]),
    id2 = as.character(df$song[j]),
    acoustic_distance = d[idx],
    species1 = as.character(df$species[i]),
    species2 = as.character(df$species[j]),
    population1 = as.character(df$population[i]),
    population2 = as.character(df$population[j]),
    individual1 = as.character(df$Individual[i]),
    individual2 = as.character(df$Individual[j]),
    geo_distance = haversine_km(df$Lat[i], df$Lon[i], df$Lat[j], df$Lon[j]),
    stringsAsFactors = FALSE
  )

  out$pop1 <- paste(out$species1, out$population1, sep = ".")
  out$pop2 <- paste(out$species2, out$population2, sep = ".")

  out
}

# list with the populations (localities) in which each species was sampled
species_occurrence <- function(df) {
  occ <- unique(data.frame(
    species = as.character(df$species),
    population = as.character(df$population)
  ))
  lapply(split(occ$population, occ$species), sort)
}

# ---- Step 1: population structure -----------------------------------------
# conspecific pairs of songs from different individuals; same_population = 1
# when both songs come from the same locality
population_structure_data <- function(pairs) {

  x <- pairs[pairs$species1 == pairs$species2 &
               pairs$individual1 != pairs$individual2, ]

  x$species <- factor(x$species1)
  x$same_population <- factor(
    as.integer(x$population1 == x$population2),
    levels = c(0, 1)
  )

  x$individual1 <- factor(x$individual1)
  x$individual2 <- factor(x$individual2)
  x$pop1 <- factor(x$pop1)
  x$pop2 <- factor(x$pop2)

  x
}

# descriptive summary of range-wide variation per species
population_structure_summary <- function(x, digits = 3) {

  do.call(rbind, lapply(split(x, x$species, drop = TRUE), function(s) {

    within  <- s$acoustic_distance[s$same_population == "1"]
    between <- s$acoustic_distance[s$same_population == "0"]

    data.frame(
      species = as.character(s$species[1]),
      populations = length(unique(c(s$population1, s$population2))),
      individuals = length(unique(c(as.character(s$individual1), as.character(s$individual2)))),
      within_population_pairs = length(within),
      between_population_pairs = length(between),
      mean_within = round(mean(within), digits),
      mean_between = if (length(between) > 0) round(mean(between), digits) else NA,
      between_within_ratio = if (length(between) > 0) round(mean(between) / mean(within), digits) else NA,
      max_geo_distance_km = round(max(s$geo_distance), 0)
    )
  }))
}

# ---- Step 3: population shifts ----------------------------------------------
# focal-congener combinations in which the focal species has at least one
# population that co-occurs with the congener and at least one that does not
# (the only combinations in which a population shift can be tested)
shift_combinations <- function(occ) {

  sp <- names(occ)
  combos <- expand.grid(focal = sp, congener = sp, stringsAsFactors = FALSE)
  combos <- combos[combos$focal != combos$congener, ]

  info <- lapply(seq_len(nrow(combos)), function(k) {
    fp <- occ[[combos$focal[k]]]
    cp <- occ[[combos$congener[k]]]
    co <- fp %in% cp
    # co-occurring focal populations that can be compared with congener
    # songs recorded at another locality (same-locality pairs are excluded)
    co_ref <- vapply(fp[co], function(p) any(cp != p), logical(1))
    data.frame(
      focal = combos$focal[k],
      congener = combos$congener[k],
      cooccurring_populations = paste(fp[co], collapse = ", "),
      non_cooccurring_populations = paste(fp[!co], collapse = ", "),
      congener_populations = paste(cp, collapse = ", "),
      testable = any(co) && any(!co) && any(co_ref)
    )
  })

  do.call(rbind, info)
}

# heterospecific focal x congener pairs, oriented so that side 1 is always
# the focal species. Pairs recorded at the same locality (the sympatric
# comparisons already used in Step 2) are excluded, so the congener songs act
# as an external reference: the question is whether focal populations that
# live with the congener are shifted towards / away from the congener's song
# when compared with congener songs recorded elsewhere
population_shift_data <- function(pairs, focal, congener, occ) {

  x <- pairs[(pairs$species1 == focal & pairs$species2 == congener) |
               (pairs$species1 == congener & pairs$species2 == focal), ]

  swap <- x$species1 != focal
  for (v in c("id", "species", "population", "individual", "pop")) {
    a <- paste0(v, "1")
    b <- paste0(v, "2")
    tmp <- x[swap, a]
    x[swap, a] <- x[swap, b]
    x[swap, b] <- tmp
  }

  x <- x[x$population1 != x$population2, ]

  x$cooccur <- factor(
    as.integer(x$population1 %in% occ[[congener]]),
    levels = c(0, 1)
  )

  x$focal_population <- x$population1
  x$acoustic_distance_sc <- as.numeric(scale(x$acoustic_distance))

  x$individual1 <- factor(x$individual1)
  x$individual2 <- factor(x$individual2)
  x$pop1 <- factor(x$pop1)
  x$pop2 <- factor(x$pop2)

  x
}

# co-occurrence effect expressed relative to the typical difference among
# focal populations (posterior of b_cooccur1 / sd_pop1__Intercept)
shift_relative_to_population_sd <- function(fit, digits = 3) {

  dr <- as_draws_df(fit)
  ratio <- dr$b_cooccur1 / dr$sd_pop1__Intercept

  data.frame(
    estimate = round(median(ratio), digits),
    `l-95% UI` = round(quantile(ratio, 0.025), digits),
    `u-95% UI` = round(quantile(ratio, 0.975), digits),
    check.names = FALSE
  )
}

# regularizing priors for the within-species models (Steps 1 and 3), the
# same set used for the trait-level sympatry models: the response is
# standardized (mean 0, SD 1), so a coefficient near |1| would already be an
# implausibly large effect. The sd prior applies to every random effect,
# including the multi-membership terms
priors_ws <- c(
  prior(normal(0, 1),   class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(exponential(1), class = "sigma")
)

# prior predictive check plot (same layout as the sympatry-model checks)
prior_check_plot <- function(fit, title) {
  pp_check(fit, type = "dens_overlay", ndraws = 50) +
    coord_cartesian(xlim = c(-10, 10)) +
    labs(
      title = title,
      subtitle = "y: observed standardized acoustic distance; y_rep: draws from the priors only"
    )
}
