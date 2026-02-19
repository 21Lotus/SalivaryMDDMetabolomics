suppressPackageStartupMessages({
  library(circlize)
  library(data.table)
  library(dplyr)
})

# vectorized deterministic hash -> [0,1]
hash01 <- function(x) {
  x <- as.character(x)
  vapply(x, function(s) {
    raw <- charToRaw(s)
    sm <- sum(as.integer(raw))
    ((sm %% 9973) / 9973)
  }, numeric(1))
}

# draw a clean vertical colorbar in the right legend column
draw_colorbar <- function(col_fun, lo, mid, hi,
                          x_left, y_bot, w, h, n = 260,
                          title = "Dz") {
  ys <- seq(y_bot, y_bot + h, length.out = n + 1)
  vals <- seq(lo, hi, length.out = n)
  
  for (i in seq_len(n)) {
    rect(x_left, ys[i], x_left + w, ys[i + 1], col = col_fun(vals[i]), border = NA)
  }
  # outline whole bar
  rect(x_left, y_bot, x_left + w, y_bot + h, border = "black", lwd = 0.7)
  
  # black boxed band at 0 (white/no-change visibility)
  if (mid >= lo && mid <= hi && is.finite(lo) && is.finite(hi) && lo != hi) {
    mid_y <- y_bot + h * ((mid - lo) / (hi - lo))
    band <- h * 0.018
    rect(x_left, mid_y - band, x_left + w, mid_y + band, border = "black", lwd = 1.4)
  }
  
  # ticks
  segments(x_left + w, y_bot, x_left + w + w*0.25, y_bot, lwd = 0.8)
  segments(x_left + w, y_bot + h/2, x_left + w + w*0.25, y_bot + h/2, lwd = 0.8)
  segments(x_left + w, y_bot + h, x_left + w + w*0.25, y_bot + h, lwd = 0.8)
  
  text(x_left + w + w*0.35, y_bot, sprintf("%.3f", lo), adj = 0, cex = 0.8)
  text(x_left + w + w*0.35, y_bot + h/2, sprintf("%.3f", mid), adj = 0, cex = 0.8)
  text(x_left + w + w*0.35, y_bot + h, sprintf("%.3f", hi), adj = 0, cex = 0.8)
  
  text(x_left, y_bot + h + h*0.08, title, adj = 0, cex = 0.95)
}

make_plot <- function() {
  
  #compact page: reserve a right column for legends, avoid white space
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  plot.new()
  rect(0, 0, 1, 1, col = "#F6F7F9", border = NA)  # subtle background
  
  #read inputs
  B      <- fread("intensity_by_group.csv")
  mcat   <- fread("metabolite_category.csv")
  map    <- fread("hmdb_to_pathway.csv")
  chords <- fread("chords.csv")
  legend_df <- fread("roman_legend.csv")
  
  setnames(B, names(B), trimws(names(B)))
  setnames(mcat, names(mcat), trimws(names(mcat)))
  setnames(map, names(map), trimws(names(map)))
  setnames(chords, names(chords), trimws(names(chords)))
  setnames(legend_df, names(legend_df), trimws(names(legend_df)))
  
  if (!("hmdb_id" %in% names(B)) && ("metabolite" %in% names(B))) setnames(B, "metabolite", "hmdb_id")
  if (!("hmdb_id" %in% names(B))) stop("intensity_by_group.csv must have hmdb_id or metabolite column.")
  
  if (!("delta" %in% names(B))) {
    if (all(c("Control","Depressed") %in% names(B))) {
      B[, delta := Depressed - Control]
    } else stop("Need delta OR both Control and Depressed.")
  }
  
  df_m <- B[, .(hmdb_id, value = delta)] %>%
    left_join(mcat, by = c("hmdb_id"="hmdb_id")) %>%
    mutate(class = ifelse(is.na(class) | class=="", "Unknown", class)) %>%
    arrange(class, hmdb_id)
  
  # ---- robust scaling to show small changes ----
  vals <- df_m$value
  lo <- as.numeric(quantile(vals, 0.05, na.rm = TRUE))
  hi <- as.numeric(quantile(vals, 0.95, na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
    lo <- min(vals, na.rm = TRUE)
    hi <- max(vals, na.rm = TRUE)
  }
  mid <- 0
  col_fun <- colorRamp2(c(lo, mid, hi), c("blue", "white", "red"))
  
  #MET sectors by class (width = count)
  met_sector_df <- df_m %>%
    group_by(class) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(sector = paste0("MET_", class), xstart = 0, xend = n)
  
  #PW sectors: roman numerals, width proportional to mapped count
  counts_pw <- map %>%
    group_by(pathway_id) %>%
    summarise(n = n(), .groups = "drop") %>%
    right_join(legend_df, by = c("pathway_id"="pathway_id")) %>%
    mutate(n = ifelse(is.na(n), 1, n))
  
  pw_sector_df <- counts_pw %>%
    mutate(sector = paste0("PW_", pathway_id), xstart = 0, xend = n)
  
  #interleave MET and PW to spread numerals around circle
  met_list <- met_sector_df %>% arrange(desc(n)) %>% as.data.frame()
  pw_list  <- pw_sector_df  %>% arrange(desc(xend)) %>% as.data.frame()
  
  sector_df <- data.frame(sector=character(), xstart=numeric(), xend=numeric(), stringsAsFactors=FALSE)
  i <- 1; j <- 1
  while (i <= nrow(met_list) || j <= nrow(pw_list)) {
    if (i <= nrow(met_list)) { sector_df <- rbind(sector_df, met_list[i, c("sector","xstart","xend")]); i <- i+1 }
    if (j <= nrow(pw_list))  { sector_df <- rbind(sector_df, pw_list[j,  c("sector","xstart","xend")]); j <- j+1 }
  }
  
  #circos: reserve right column by shrinking the circos canvas
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = 1.0,
    cell.padding = c(0.0005, 0, 0.0005, 0),
    track.margin = c(0.0005, 0.0005),
    canvas.xlim = c(-0.1, 0.70),   # <- circos occupies left, right is legend column
    canvas.ylim = c(-1.12, 1.12)
  )
  circos.initialize(factors = sector_df$sector, xlim = sector_df[, c("xstart","xend")])
  
  #PW labels ONLY (roman numerals)
  circos.trackPlotRegion(ylim=c(0,1), track.height=0.06, bg.border=NA,
                         panel.fun=function(x,y){
                           s <- get.cell.meta.data("sector.index")
                           if (!startsWith(s,"PW_")) return()
                           lab <- sub("^PW_","",s)
                           circos.text(mean(get.cell.meta.data("xlim")), 0.5, lab,
                                       cex=0.58, facing="clockwise", niceFacing=TRUE)
                         }
  )
  
  #heatmap ring (thin)
  circos.trackPlotRegion(ylim=c(0,1), track.height=0.045, bg.border=NA,
                         panel.fun=function(x,y){
                           s <- get.cell.meta.data("sector.index")
                           if (!startsWith(s,"MET_")) return()
                           cls <- sub("^MET_","",s)
                           mets <- df_m$hmdb_id[df_m$class==cls]
                           v <- df_m$value[match(mets, df_m$hmdb_id)]
                           for (k in seq_along(mets)) {
                             circos.rect(k-1, 0, k, 1, col=col_fun(v[k]), border=NA)
                           }
                         }
  )
  
  #class ring (thin)
  classes <- met_sector_df$class
  pal <- grDevices::hcl.colors(length(classes), "Set 3")
  class_cols <- setNames(pal, classes)
  
  circos.trackPlotRegion(ylim=c(0,1), track.height=0.02, bg.border=NA,
                         panel.fun=function(x,y){
                           s <- get.cell.meta.data("sector.index")
                           if (!startsWith(s,"MET_")) return()
                           cls <- sub("^MET_","",s)
                           xl <- get.cell.meta.data("xlim")
                           circos.rect(xl[1],0,xl[2],1,col=class_cols[cls], border=NA)
                         }
  )
  
  #thicker, denser, across-circle links
  met_to_sector <- df_m %>%
    group_by(class) %>%
    mutate(pos = row_number()) %>%
    ungroup() %>%
    transmute(hmdb_id, sector1=paste0("MET_",class), x1=pos-0.5, class=class)
  
  pw_width_map <- setNames(pw_sector_df$xend, pw_sector_df$pathway_id)
  
  link_df <- chords %>%
    rename(hmdb_id=from, pathway_id=to) %>%
    left_join(met_to_sector, by="hmdb_id") %>%
    mutate(
      sector2 = paste0("PW_", pathway_id),
      pw_n = pw_width_map[pathway_id],
      # spread landing positions across entire PW arc
      x2 = pmax(0.6, pmin(as.numeric(pw_n) - 0.6, hash01(hmdb_id) * as.numeric(pw_n)))
    ) %>%
    filter(!is.na(sector1), !is.na(sector2), is.finite(x2))
  
  # draw links (dense look)
  apply(link_df, 1, function(r){
    cls <- r[["class"]]
    circos.link(
      r[["sector1"]], as.numeric(r[["x1"]]),
      r[["sector2"]], as.numeric(r[["x2"]]),
      col = adjustcolor(class_cols[cls], alpha.f = 0.60),
      border = NA,
      lwd = 2.2,
      h.ratio = 0.1   # deep arcs across circle
    )
  })
  
  #clean legends
  par(xpd = NA)
  
  #convert device coordinates (0–1) to user coordinates
  to_user_x <- function(x) grconvertX(x, from="ndc", to="user")
  to_user_y <- function(y) grconvertY(y, from="ndc", to="user")
  
  #roman numeral legend
  legend(
    x = to_user_x(0.69),
    y = to_user_y(0.76),
    legend = paste0(legend_df$pathway_id, "  ", legend_df$pathway_name),
    bty = "n",
    cex = 0.85,
    xjust = 0,
    yjust = 1
  )
  
  #Dz color bar
  #position box in device space
  x1 <- to_user_x(0.73)
  x2 <- to_user_x(0.77)
  y1 <- to_user_y(0.24)
  y2 <- to_user_y(0.44)
  
  draw_colorbar(
    col_fun,
    lo,
    0,
    hi,
    x_left = x1,
    y_bot = y1,
    w = (x2 - x1),
    h = (y2 - y1),
    n = 260,
    title = "Dz"
  )
  
  par(xpd = FALSE)
}

#render in pane
make_plot()

#export (tight)
pdf("circos_landscape.pdf", width = 20, height = 6.8)
make_plot()
dev.off()

png("circos_landscape.png", width = 4200, height = 2400, res = 300)
make_plot()
dev.off()

cat("Saved: circos_landscape.pdf and circos_landscape.png\n")