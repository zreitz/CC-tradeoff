library(jsonlite)

writeConfig <- function(path, desc, 
                        replicates = 1, 
                        rounds = 100000, 
                        write_every = 10, 
                        seed = 93117,
                        pop_size = 10000,
                        patch_lifetime = 10000,
                        mort_to_disp = 1.21836,
                        uniform_fraction = 0,
                        max_biomass = 121.012,
                        mut_prob = 0,
                        mut_stdev = 0.001,
                        mut_range = c(0, 0.924148),
                        initial = list(
                          list(alpha = 0.211, patches = 1000)
                        )) {
  
  # Construct the list structure
  config <- list(
    description = desc,
    replicates = replicates,
    rounds = rounds,
    write_every = write_every,
    seed = seed,
    pop_size = pop_size,
    patch_lifetime = patch_lifetime,
    mort_to_disp = mort_to_disp,
    uniform_fraction = uniform_fraction,
    max_biomass = max_biomass,
    mut_prob = mut_prob,
    mut_stdev = mut_stdev,
    mut_range = mut_range,
    initial = initial
  )
  
  # Write to file
  # auto_unbox = TRUE: prevents single values from being wrapped in [brackets]
  # pretty = TRUE: makes the JSON human-readable
  write_json(config, path, auto_unbox = TRUE, pretty = TRUE)
}

comparePlot <- function(points, lims = c(0, 50000)) {
  points <- points[points$time > lims[1] * 0.9 & points$time < lims[2] * 1.1, ]
  ggplot(points, aes(x = time, y = alpha)) +
    geom_hline(yintercept = 0.0738069, linetype = "dashed") +
    geom_point(size = 0.4, show.legend = F, color = "black") +
    # geom_hline(yintercept = 0.0743315, linetype = "dashed") +
    # geom_hline(yintercept = 0.226797, linetype = "dashed") +
    labs(
      x = NULL,
      y = NULL,
      color = "Patch\nOccupancy"
    ) +
    # scale_color_gradientn(colors = tol(23)[10:23],
    #                       trans = "log",
    #                       limits = c(0.0001, 1),
    #                       breaks = c(0.0001, 0.001, 0.01, 0.1, 1),
    #                       labels = c("0.0001", 0.001, 0.01, 0.1, 1)) +
    #coord_cartesian(xlim = lims, ylim = c(0.05, 0.6)) +
    scale_x_continuous(
      limits = c(lims[1] - 1000, lims[2] + 1000),
      expand = c(0,0),
      breaks = lims,
      labels = sprintf("%sK", round(lims/1000, 1))
      ) +
    scale_y_continuous(
      limits = c(0.05, 0.61),
      expand = c(0,0),
      breaks = c(0.1, 0.3, 0.5),
      labels = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      aspect.ratio = 1.7
    )
}