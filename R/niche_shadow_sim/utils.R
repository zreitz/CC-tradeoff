
writeConfig <- function(
  path,
  desc,

  mortality = 1.5,
  mut_stdev = 5e-4,
  mut_uniform = 0.05,

  write_every = 10,
  rounds = 1000000,
  initial = 0.227
) {
  txt <- sprintf('{
  "description": "%s",
  "rounds": %i,
  "write_every": %i,
  "mortality": %s,
  "max_biomass": 5.35,
  "mut_stdev": %s,
  "mut_uniform": %s,
  "mut_range": [
    0.0738069,
    0.822162
  ],
  "seed": 93117,
  "initial": [
    %s
  ]
}', desc, rounds, write_every, mortality, mut_stdev, mut_uniform, initial)
  write(txt, file = path)
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