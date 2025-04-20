# Real-life df implementation 
require(effectFusion)

df <- read.csv("Clean_Dataset.csv")
df <- df[,-1]

# Create combined journey and time columns
df$journey <- paste(df$source_city, df$destination_city, sep = "-")
df$time <- paste(df$departure_time, df$arrival_time, sep = "-")

# Divide into bins of 5 hrs
df$duration_bin <- cut(
  df$duration,
  breaks = seq(0, max(df$duration, na.rm = TRUE) + 5, by = 5),
  right = FALSE,
  include.lowest = TRUE
)
df$days_left_bin <- cut(
  df$days_left,
  breaks = seq(0, max(df$days_left, na.rm = TRUE) + 5, by = 5),
  right = FALSE,
  include.lowest = TRUE
)

# Select relevant predictors
X <- df[c("airline", "journey", "stops", "class", "duration_bin", "days_left_bin")]
y <- df$price

# Encode variables properly
X$airline         <- factor(X$airline)                                   # nominal
X$journey         <- factor(X$journey)                                   # nominal
X$stops           <- factor(X$stops, levels = c("zero", "one", "two_or_more"), ordered = TRUE)
X$class           <- factor(X$class, levels = c("Economy", "Business"), ordered = TRUE)
X$duration_bin    <- factor(X$duration_bin, ordered = TRUE)
X$days_left_bin   <- factor(X$days_left_bin, ordered = TRUE)

# Scale target
y <- scale(y)

# Define type list: "n" = nominal, "o" = ordinal
type_list <- c("n", "n", "o", "o", "o", "o")

# Fit model
fusion_fit <- effectFusion(
  y = y,
  X = X,
  types = type_list,
  method = "SpikeSlab",
  prior = list(r = 20000, G0 = 20),
  mcmc = list(M = 10000, burnin = 5000),
  modelSelection = "binder"
)

invisible_output <- capture.output(coef_det <- summary(fusion_fit)[2:nrow(summary(fusion_fit)),])
fit_coef <- coef_det[,1]
summary <- summary(fusion_fit)
print(summary)


save(fit_coef, summary, file = "Real-Life-Study.RData")
