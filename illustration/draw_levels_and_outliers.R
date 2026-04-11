library(ggplot2)
library(dplyr)
library(scales)
library(grid)

# ------------------------------------------------------------
# helper: build plotting data
# ------------------------------------------------------------
build_flagged_plot_data <- function(data, varname, idx_huber, n_segments = 4){
  tt = index(data)
  yy = as.numeric(data[, varname])
  
  df = data.frame(
    time = as.Date(tt),
    value = yy,
    is_flag = idx_huber == 0
  )
  
  n = nrow(df)
  break_pts = round(seq(1, n + 1, length.out = n_segments + 1))
  
  df$segment = NA_integer_
  for(s in seq_len(n_segments)){
    i_start = break_pts[s]
    i_end   = break_pts[s + 1] - 1
    df$segment[i_start:i_end] = s
  }
  
  df$segment = factor(df$segment, levels = 1:n_segments)
  df
}

# ------------------------------------------------------------
# helper: four interior date breaks for each facet
# divide the span into 5 equal parts and take 4 interior points
# ------------------------------------------------------------
facet_date_breaks <- function(x){
  x0 = min(x)
  x1 = max(x)
  
  if(x0 == x1){
    return(x0)
  }
  
  x_num = as.numeric(c(x0, x1))
  b = seq(x_num[1], x_num[2], length.out = 6)
  as.Date(round(b[2:5]), origin = "1970-01-01")
}

# ------------------------------------------------------------
# plotting function
# ------------------------------------------------------------
plot_flagged_series <- function(data, varname, idx_huber, ylab_name){
  
  df = build_flagged_plot_data(data, varname, idx_huber, 4)
  df_flag = df %>% filter(is_flag)
  
  ggplot(df, aes(x = time, y = value)) +
    
    geom_line(linewidth = 0.35, colour = "black") +
    
    geom_vline(
      data = df_flag,
      aes(xintercept = time),
      linewidth = 0.25,
      colour = "grey50"
    ) +
    
    geom_point(
      data = df_flag,
      size = 0.7,
      colour = "black"
    ) +
    
    facet_wrap(~ segment, ncol = 1, scales = "free_x") +
    
    scale_x_date(
      breaks = facet_date_breaks,
      date_labels = "%Y-%m",
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    labs(
      x = NULL,
      y = ylab_name
    ) +
    
    theme_bw(base_size = 10) +
    theme(
      strip.text = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey80"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(0.2, "lines"),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10),
      plot.margin = margin(4, 6, 4, 6)
    )
}

# ------------------------------------------------------------
# generate plots
# ------------------------------------------------------------
p1 = plot_flagged_series(
  data,
  varnames[1],
  idx_huber,
  "Jökulsá"
)

p2 = plot_flagged_series(
  data,
  varnames[2],
  idx_huber,
  "Vatndalsá"
)

p1
p2

# ------------------------------------------------------------
# save figures as PNG
# ------------------------------------------------------------

ggsave(
  filename = "jokulsa_flags.png",
  plot = p1,
  width = 6,
  height = 4,
  dpi = 300
)

ggsave(
  filename = "vatndalsa_flags.png",
  plot = p2,
  width = 6,
  height = 4,
  dpi = 300
)
