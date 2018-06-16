# ----------------------------------------------------------------------------
# RLog - Tactrix Standalone Logger Analysis Tool for the DIY NC-MX5 Community
# ----------------------------------------------------------------------------
# This utility is inspired in it's entirety by the ideas / efforts of 
# Skyavonne on miata.net's forum. It's a very interesting way to think
# about tuning, and gave me a opporunity to explore R quite a bit :)
#
#
#
# ----------------------------------------------------------------------------
# HOW TO RUN RLOG SCRIPT
# ----------------------------------------------------------------------------
# OS-X: Command + Option + R
# Win: Control + Alt + R 
#
# Notes:
# RLog will attempt install two packages if not already installed upon
# execution. The "Purrr" library will specifically ask if you want to install 
# using sources that need compilation. Just enter "n" when prompted, it will 
# save time. Once the packages are installed, you won't be prompted again.
#
#
#
# ----------------------------------------------------------------------------
# [!] BELOW IS THE ONLY LINE THAT REQUIRES CHANGE
# ----------------------------------------------------------------------------
# log_path can refer to a directory of logs or an individual log file. 
# Logs must be csv formated and should adhere to supplied logcfg.txt

log_path = "[!] REPLACE WITH FILE PATH TO LOG DIR OR LOG FILE [!]"




# REQUIRED PACKAGES

if(!require(plyr)){
  install.packages("plyr")
  library(plyr)
}

if(!require(purrr)){
  install.packages("purrr")
  library(purrr)
}




# HELPER FUNCTIONS

# math helper functions
math.midpoint <- function(x, y) x + ((y - x) * 0.5)
math.round.2  <- function(x) round(x,2)
math.round.3  <- function(x) round(x,3)
math.round.4  <- function(x) round(x,4)

util.count <- function(df) length(df[is.finite(df)])
util.count.above.threshold <- function(df, thresh) length(df[df>thresh])
util.count.increasing <- function(a) {
  cur_count  <- 0
  prev_value <- 0
  for (i in seq_along(a)) {
    if (is.finite(a[i]) & a[i] > prev_value) (cur_count <- cur_count + 1)
    prev_value <- a[i]
  }
  cur_count
}




# FILTERS

# base filters
fitler.ol.all        <- function(df) subset(df, FUELSYS == 4 | FUELSYS == 1)
filter.ol.trans      <- function(df) subset(df, FUELSYS == 1)
filter.ol            <- function(df) subset(df, FUELSYS == 4)
filter.cl            <- function(df) subset(df, FUELSYS == 2 | FUELSYS == 16)
filter.ect           <- function(df, from = -100, to = 200) subset(df, from <= ECT & ECT <= to)
filter.iat           <- function(df, from = -100, to = 200) subset(df, from <= IAT & IAT <= to)
filter.app           <- function(df, from = 0, to = 100) subset(df, from <= APP & APP <= to)
filter.vss           <- function(df, from = 0, to = 300) subset(df, from <= VSS & VSS <= to)
filter.ipw           <- function(df, from = 0, to = 50) subset(df, from <= FUEL_PW & FUEL_PW <= to)
filter.kr            <- function(df, from = 0, to = 10) subset(df, from <= KNOCKR & KNOCKR <= to)
filter.egr           <- function(df, from = 0, to = 52) subset(df, form <= EGRP_STEPS & EGRP_STEPS <= to)
filter.etc           <- function(df, from = 0, to = 100) subset(df, from <= ETC_DSD & ETC_DSD <= to)

# composable filters
filter.app.off       <- function(df) filter.app(df, to = 0)
filter.app.on        <- function(df) filter.app(df, from = 2.3)
filter.app.wot       <- function(df) filter.app(df, from = 82)
filter.ect.warm      <- function(df) filter.ect(df, from = 80)
fitler.ect.cold      <- function(df) filter.ect(df, to = 80)
filter.injectors.off <- function(df) filter.ipw(df, from = 0, to = 0)
filter.injectors.on  <- function(df) filter.ipw(df, from = 0.001)
filter.vss.stop      <- function(df) filter.vss(df, to = 0)




# BINNING AND BREAKPOINTS

# generic binning function
bin.data <- function(data, field, breaks) cut(data[[field]], breaks = breaks, dig.lab = 10)

# function creates bins that evenly surround a desired mid point
bin.ranges.create <- function(from, to, by) (seq(from - (by * 0.5), to + (by * 0.5), by = by))

# converts binned ranges back to midpoints
bin.ranges.to.midpoints <- function(bins) {
  from <- as.numeric(sub("\\((.+),.*", "\\1", bins))
  to   <- as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", bins))
  midpoints <- mapply(math.midpoint, from, to)
  midpoints <- sapply(midpoints, math.round.3)
  midpoints
}

# breakpoints used for binning and plot points
breakpoints.APP          <- c(0,2.3,5,7.5,10,15,20,25,30,35,40,45,50,60,70,80,90,100)
breakpoints.ECT          <- bin.ranges.create(from = -40, to = 120, by = 10)
breakpoints.EGRP_STEPS   <- bin.ranges.create(from = 0, to = 50, by = 5)
breakpoints.EQ_RATIO_ACT <- bin.ranges.create(from = 0.7, to = 1.4, by = 0.1)
breakpoints.ETC_DSD      <- bin.ranges.create(from = 0, to = 100, by = 12.5)
breakpoints.LOAD         <- bin.ranges.create(from = 0, to = 1, by = 0.0625)
breakpoints.IAT          <- bin.ranges.create(from = -40, to = 120, by = 10)
breakpoints.IDC          <- bin.ranges.create(from = 0, to = 100, by = 10)
breakpoints.MAF          <- c(0,15,30,200)
breakpoints.RPM          <- bin.ranges.create(from = 0, to = 7500, by = 500)
breakpoints.SPARKADV     <- bin.ranges.create(from = -20, to = 50, by = 5)
breakpoints.VPWR         <- bin.ranges.create(from = 5.25, to = 16.50, by = 1.25)

# returns a vector with breakpoints based on midpoints
# uneven intervals indicate breakpoints not based on seq function
center.breakpoints <- function(breakpoints) {
  h = head(breakpoints,1)
  t = tail(breakpoints,1)
  l = length(breakpoints)-1
  i = abs((t - h) / l)
  s = h+(i*0.5)
  e = t-(i*0.5)
  test_breakpoints <- bin.ranges.create(from = s, to = e, by = i)
  if (all(breakpoints == test_breakpoints))
    return (seq(from = s, to = e, by = i))
  return (breakpoints)
}  

# specific PIDs we want binned for plots and matrices.
# These can be expended to include new PIDS, but please
# note the format "bin.PID." The PID must match 
# associated column name from CSV.
bin.APP         <- function(data) bin.data(data, 'APP'         , breakpoints.APP)
bin.ECT         <- function(data) bin.data(data, 'ECT'         , breakpoints.ECT)
bin.EGRP_STEPS  <- function(data) bin.data(data, 'EGRP_STEPS'  , breakpoints.EGRP_STEPS)
bin.EQ_RATIO_ACT<- function(data) bin.data(data, 'EQ_RATIO_ACT', breakpoints.EQ_RATIO_ACT)
bin.ETC_DSD     <- function(data) bin.data(data, 'ETC_DSD'     , breakpoints.ETC_DSD)
bin.LOAD        <- function(data) bin.data(data, 'LOAD'        , breakpoints.LOAD)
bin.IAT         <- function(data) bin.data(data, 'IAT'         , breakpoints.IAT)
bin.IDC         <- function(data) bin.data(data, 'IDC'         , breakpoints.IDC)
bin.MAF         <- function(data) bin.data(data, 'MAF'         , breakpoints.MAF)
bin.RPM         <- function(data) bin.data(data, 'RPM'         , breakpoints.RPM)
bin.SPARKADV    <- function(data) bin.data(data, 'SPARKADV'    , breakpoints.SPARKADV)
bin.VPWR        <- function(data) bin.data(data, 'VPWR'        , breakpoints.VPWR)

# utility that converts a string into a bin functions
# e.g, "APP" as an input would evaluate to the bin.APP function
string.to.bin.func <- function(x) {
  x <- paste("bin", x, sep=".")
  x <- eval(parse(text = x))
}

# utility that converts a string into a breakpoints variable
string.to.breakpoints <- function(x) {
  x <- paste("breakpoints", x, sep=".")
  x <- eval(parse(text = x))
}




# MATRIX FUNCTIONS

matrix.divide <- function(a,b){
  m <- a/b
  m[is.nan(m)] <- NA
  m <-apply(m, 2, math.round.3)
  m
}

# generic matrix creation function
matrix.create <- function(data, x, y, f) {
  bin.col <- string.to.bin.func(x)
  bin.row <- string.to.bin.func(y)
  df          <- ddply(data, .(bin.col(data), bin.row(data)), .drop = F, .fun = f)
  df$V1       <- sapply(df$V1, math.round.3)
  col_labels  <- levels(df$`bin.col`)
  row_labels  <- levels(df$`bin.row`)
  m           <- matrix(df$V1, length(row_labels), length(col_labels), byrow = F)
  rownames(m) <- row_labels
  colnames(m) <- col_labels
  m[!is.finite(m)] <- NA
  m
}

# generic matrix creation function but using midpoints for the axes
matrix.create.with.midpoints <- function(data, x, y, f) {
  m <- matrix.create(data, x, y, f)
  rownames(m) <- bin.ranges.to.midpoints(rownames(m))
  colnames(m) <- bin.ranges.to.midpoints(colnames(m))
  m
}

matrix.create.with.row.midpoints <- function(data, x, y, f) {
  m <- matrix.create(data, x, y, f)
  rownames(m) <- bin.ranges.to.midpoints(rownames(m))
  m
}

matrix.create.with.col.midpoints <- function(data, x, y, f) {
  m <- matrix.create(data, x, y, f)
  colnames(m) <- bin.ranges.to.midpoints(colnames(m))
  m
}




# PLOTTING FUNCTIONS

scatterplot.xy <- function(df, x, y, col = '#3d6870', ...) {
  xrange <- center.breakpoints(string.to.breakpoints(x))
  yrange <- center.breakpoints(string.to.breakpoints(y))
  par(bg = '#333333')
  par(fg = '#777777')
  plot(df[[x]],
       df[[y]], 
       col.main = '#cccccc',
       xlab     = x, 
       ylab     = y,
       xlim     = c(head(xrange,1),tail(xrange,1)),
       ylim     = c(head(yrange,1),tail(yrange,1)),
       col      = col,
       cex      = 0.5,
       col.axis = '#aaaaaa',
       col.lab  = '#aaaaaa',
       cex.axis = 0.9,
       cex.lab  = 0.9,
       mgp      = c(2,.5,0), ...)
  grid (NULL, NULL, col = '#777777', lwd = 0.5, lty = 3)
}

scatterplot.overlay.mean <- function(data, x, y) {
  bin.col    <- string.to.bin.func(x)
  mean_df    <- ddply(data, .(bin.col(data)), .drop = F, .fun = function(d) mean(d[[y]]))
  mean_df$V1 <- sapply(mean_df$V1, math.round.3)
  mean_df$`bin.col(data)`<- bin.ranges.to.midpoints(mean_df$`bin.col(data)`)
  lines(mean_df, type="o", col = "#ffffff", lty = 1, lwd = 1)
}

scatterplot.overlay.quantile <- function(data, x ,y, q1, q2, q3) {
  bin.col <- string.to.bin.func(x)
  overlay.create <- function(q,line_col, line_type, line_width, lab_col, lab_pos) {
    pct_df    <- ddply(data, .(bin.col(data)), .drop = F, .fun = function(d) quantile(d[[y]], probs = q, names = F))
    pct_df$V1 <- sapply(pct_df$V1, math.round.3)
    pct_df$`bin.col(data)` <- bin.ranges.to.midpoints(pct_df$`bin.col(data)`)
    lines(pct_df, type='o', col = line_col, lty = line_type, lwd = line_width)
    if (length(pct_df$V1) != 0) {
      text(pct_df, labels = pct_df$V1, cex = .7, font=1, col = lab_col, pos = lab_pos)
    }
  }
  overlay.create(q1, '#ff816e', line_type = 2, line_width = 1, '#e7ba23', 1)
  overlay.create(q2, '#a3e2ff', line_type = 1, line_width = 1, '#eeeeee', 2)
  overlay.create(q3, '#ff816e', line_type = 2, line_width = 1, '#e7ba23', 3)
}

# scatter plot function with overlays for mean, median, definable percentiles)
scatterplot <- function(df, x, y, overlays = T, ...) {
  scatterplot.xy(df, x, y, ...)
  if (overlays) scatterplot.overlay.quantile(df, x, y, 0.05, 0.50, 0.95)
}




# FILE OUTPUT RELATED 
output.file.path.add <- function(filename) {
  cat("[output]:", filename,"\n")
  output_dir <- 'analysis'
  if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = T)
  file.path(output_dir,filename)
}

png.set <- function(filename = filename, ...) {
  png(filename = output.file.path.add(filename), width = 1920, height = 1080, res = 144)
}



# LOGS

log.runtime     <- function(data) tail(data$time,1)-head(data$time,1)
log.distance    <- function(data) log.runtime(data) / 60 * mean(data$VSS) / 60 * 0.62
log.samplecount <- function(data) tail(data$sample,1)-head(data$sample,1)

# takes a log and filter, returns a summary data frame
log.summary     <- function(df) {
  battv_mean    <- mean(df$VPWR) 
  baro_min      <- min(df$BARO)
  baro_max      <- max(df$BARO)
  distance      <- log.distance(df)
  ect_min       <- min(df$ECT)
  ect_max       <- max(df$ECT)
  hidet_sw      <- max(df$HIDET_SW)
  iat_min       <- min(df$IAT)
  iat_max       <- max(df$IAT)
  idc_max_pct   <- max(df$IDC)
  kr_count      <- util.count.increasing(df$KNOCKR)
  kr_max        <- max(df$KNOCKR)
  lambda_mean   <- mean(df$EQ_RATIO_ACT)  
  ltft_mean     <- mean(df$LTFT)
  lines         <- tail(df$sample,1) - head(df$sample, 1)
  rpm_max       <- max(df$RPM)
  runtime       <- log.runtime(df)
  stft_mean     <- mean(df$STFT)
  total_ft_mean <- mean(df$COMBINED_FT)
  wot_app_count <-  util.count.above.threshold(df$APP, 82)
  ol_load_count <- util.count.above.threshold(df$LOAD, 0.75)
  output         <- cbind(runtime, distance, lines, battv_mean,
                         stft_mean, ltft_mean, total_ft_mean, lambda_mean, 
                         idc_max_pct, hidet_sw, kr_count, kr_max, ect_min, 
                         ect_max, iat_min, iat_max, baro_min, baro_max,
                         rpm_max,ol_load_count,wot_app_count)
  apply(output, 2, math.round.3)
}

# takes a list of log data frames and filter, returns a data frame with summary information
log.list.summary <- function(log_list) t(as.data.frame(lapply(log_list, log.summary)))

# merge logs into one data frame
log.list.merge <- function(log_list) do.call (rbind.data.frame, log_list)

# load individual log (csv), returns a log data frame
log.read.csv <- function(filename) {
  df <- read.csv(file = filename, header = TRUE, sep = ",")
  
  # below is a hack for a tactrix bug w/ stand alone logger when logging lots of PIDs.
  # the unusually large values are not real, looks like a bit based error, so we just 
  # zero them out for the time being.
  df$O2S12_FT[df$O2S12_FT > 7] <- 0.0  
  
  # calculated PIDs appended to log dataframe
  df$COMBINED_FT <- df$LTFT + df$STFT + df$O2S12_FT
  df$MANVAC      <- (df$BARO - df$MAP)
  df$IDC         <- df$FUEL_PW * df$RPM / 120 * 0.001 * 100
  df
}

# load directory or individual logs, returns list of log data frames
log.load <- function(path) {
  filenames <- if ( file.exists(path) && !dir.exists(path) ) {
    setwd(dirname(path))
    c(basename(path))
  } else {
    setwd(path)
    list.files(pattern = ".csv")
  }
  log_list <- lapply(filenames, log.read.csv)
  names(log_list) <- filenames
  log_list
}

# this is where we add routines to include in the analysis pass.
log.analyze  <- function(path) {
  

  # LOG SUMMARY
  
  # load logs into list of data frames
  log_list <- log.load(path)
  
  # remove fuelcuts from log data
  filtered_log_list <- lapply(log_list, compose(filter.ect.warm, filter.injectors.on))
  
  summary <- as.data.frame(log.list.summary(filtered_log_list))
  write.csv (summary, file = output.file.path.add('summary.csv'))

  # merge unfiltered logs into a single data_frame
  df_merged <- log.list.merge(log_list)
  

  

  # LOG PLOTS
  
  png.set('plot CL-OL unfiltered (load x lambda).png')
  scatterplot(df_merged,
              x    = 'LOAD', 
              y    = 'EQ_RATIO_ACT',
              main = 'Filters: NONE')
  dev.off()
  
  png.set('plot CL warm non-idle (load x lambda).png')
  scatterplot(compose(filter.cl, 
                      filter.ect.warm,
                      filter.app.on,
                      filter.injectors.on)
              (df_merged),
              x    = 'LOAD',
              y    = 'EQ_RATIO_ACT',
              main = 'Filters: closed loop, engine warm, app engaged, no fuelcut')
  dev.off()
  
  png.set('plot CL warm idle (load x lambda).png')
  scatterplot(compose(filter.cl, 
                      filter.ect.warm,
                      filter.app.off,
                      filter.vss.stop,
                      filter.injectors.on)
              (df_merged),
              x    = 'LOAD',
              y    = 'EQ_RATIO_ACT',
              main = 'Filters: closed loop, engine warm, app released, injectors on, vehicle at stop')
  dev.off()
  
  
  png.set('plot OL warm decel fuelcut (load x lambda).png')
  scatterplot(compose(filter.ol,
                      filter.ect.warm,
                      filter.app.off,
                      filter.injectors.off)
              (df_merged),
              x    = 'LOAD',
              y    = 'EQ_RATIO_ACT',
              main = 'Filters: open loop, engine warm, accel released, injectors off (fuelcut)')
  dev.off()
  
  png.set('plot OL-trans warm (load x lambda).png')
  scatterplot(compose(filter.ol.trans,
                      filter.ect.warm)
              (df_merged),
              x    = 'LOAD',
              y    = 'EQ_RATIO_ACT',
              main = 'Filters: transitional open loop, engine warm')
  dev.off()
  
  png.set('plot OL warm non-idle (load x lambda).png')
  scatterplot(compose(filter.ol,
                      filter.ect.warm,
                      filter.app.on,
                      filter.injectors.on)
              (df_merged),
              x    = 'LOAD',
              y    = 'EQ_RATIO_ACT',
              main = 'Filters: open loop, engine warm, accel engaged, injectors on')
  dev.off()
  
  png.set('plot OL warm (load x idc).png')
  scatterplot(compose(filter.ol,
                         filter.ect.warm)
              (df_merged),
              x    = 'LOAD',
              y    = 'IDC',
              main = 'Filters: open loop, engine warm')
  dev.off()
  
  png.set('plot OL warm wot (rpm x lambda).png')
  scatterplot(compose(filter.ol,
                      filter.ect.warm,
                      filter.app.wot)
              (df_merged),
              x    = 'RPM',
              y    = 'EQ_RATIO_ACT',
              main = 'Filters: open loop, engine warm, WOT')
  dev.off()
  
  png.set('plot OL warm wot (rpm x idc).png')
  scatterplot(compose(filter.ol,
                         filter.ect.warm,
                         filter.app.wot)
              (df_merged),
              x    = 'RPM',
              y    = 'IDC',
              main = 'Filters: open loop, engine warm, WOT')
  dev.off()
  
  png.set('plot CL-OL warm (load x app).png')
  scatterplot(compose(filter.ect.warm,
                      filter.app.on)
              (df_merged),
              x    = 'LOAD',
              y    = 'APP',
              main = 'Filters: engine warm, APP engaged')
  dev.off()
  
  
  png.set('plot CL-OL warm (load x etc).png')
  scatterplot(compose(filter.ect.warm,
                      filter.app.on)
              (df_merged),
              x    = 'LOAD',
              y    = 'ETC_DSD',
              main = 'Filters: engine warm, APP engaged')
  dev.off()
  
  
  png.set('plot CL-OL warm (app x etc).png')
  scatterplot(compose(filter.ect.warm,
                      filter.app.on)
              (df_merged),
              x    = 'APP',
              y    = 'ETC_DSD',
              main = 'ETC_DSD vs APP | Filters: engine warm, APP engaged')
  dev.off()
  
  
  
  
  # CLOSED LOOP MATRICES
  
  df <- compose(filter.cl, filter.ect.warm, filter.injectors.on)(df_merged)
  
  # sample count
  m_counts <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) util.count(df$sample))
  write.csv (m_counts, file = output.file.path.add('CL sample count.csv'))
    
  # SPARK
  # spark advance mean
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) mean(df$SPARKADV))
  write.csv (m, file = output.file.path.add('CL sparkadv mean (load x rpm).csv'))
    
  # knock retard max
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) suppressWarnings(max(df$KNOCKR)))
  write.csv (m, file = output.file.path.add('CL knockr max (load x rpm).csv'))
    
  # knock retard count
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) util.count.increasing(df$KNOCKR))
  write.csv (m, file = output.file.path.add('CL knockr count (load x rpm).csv'))

  # knock retard percentage (kr_count/sample_count)
  m <- matrix.divide(m, m_counts)
  write.csv (m, file = output.file.path.add('CL knockr percent (load x rpm).csv'))

  # egr steps mean
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) mean(df$EGRP_STEPS))
  write.csv (m, file = output.file.path.add('CL egr steps (load x rpm).csv'))
  
  # FUEL
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) (mean(df$COMBINED_FT) * 0.01) - (mean(summary$total_ft_mean) * 0.01))
  write.csv (m, file = output.file.path.add('CL fuel trim (load x rpm).csv'))
  
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) mean(df$EQ_RATIO_ACT))
  write.csv (m, file = output.file.path.add('CL lambda mean (load x rpm).csv'))
  
  m <- matrix.create(df,'MAF','RPM',function(df) mean(df$COMBINED_FT))
  write.csv(m, file = output.file.path.add(('CL combined FT mean (maf x rpm).csv')))
  
  
  
  
  # OPEN LOOP MATRICES
  
  df <- compose(filter.ol, filter.ect.warm)(df_merged)
  
  # sample count
  m_counts <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) util.count(df$sample))
  write.csv (m_counts, file = output.file.path.add('OL sample count.csv'))

  # SPARK
  # spark advance mean
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) mean(df$SPARKADV))
  write.csv (m, file = output.file.path.add('OL sparkadv mean (load vs rpm).csv'))
  
  # knock retard max
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) suppressWarnings(max(df$KNOCKR)))
  write.csv (m, file = output.file.path.add('OL knockr max (load vs rpm).csv'))
  
  # knock retard count
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) util.count.increasing(df$KNOCKR))
  write.csv (m, file = output.file.path.add('OL knockr count (load vs rpm).csv'))
  
  # knock retard percentage (kr_count/sample_count)
  m <- matrix.divide(m, m_counts)
  write.csv (m, file = output.file.path.add('OL knockr percent (load vs rpm).csv'))
  
  # FUEL
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) (mean(df$COMBINED_FT) * 0.01) - (mean(summary$total_ft_mean) * 0.01))
  write.csv (m, file = output.file.path.add('OL fuel trim (load vs rpm).csv'))
  
  m <- matrix.create.with.midpoints(df, 'LOAD', 'RPM', function(df) mean(df$EQ_RATIO_ACT))
  write.csv (m, file = output.file.path.add('OL lambda mean (load vs. rpm).csv'))
  
  m <- matrix.create.with.midpoints(df, 'ETC_DSD', 'RPM', function(df) mean(df$EQ_RATIO_ACT))
  write.csv (m, file = output.file.path.add('OL lambda mean (throttle vs. rpm).csv'))
  
  m <- matrix.create.with.row.midpoints(df,'MAF','RPM',function(df) mean(df$COMBINED_FT))
  write.csv(m, file = output.file.path.add(('OL combined FT mean (maf vs rpm).csv')))
  
  cat('\n>',nrow(df_merged),'log lines.')
}

log.analyze(log_path)
