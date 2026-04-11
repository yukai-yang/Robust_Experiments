# cannot be run directly
# should be sourced by other script

# read data
river <- read.table(
  "rivers.dat",
  header = TRUE,
  quote = "\"",
  stringsAsFactors = FALSE
)

# drop index column if needed
if (!"Date" %in% colnames(river)[1]) {
  river <- river[, -1]
}

# convert date
river$Date <- as.Date(river$Date, "%d-%b-%y")

# levels (log series already provided)
joku <- river$LogJoku
vatn <- river$LogVatn

# build time-indexed object
data <- zoo::zoo(cbind(joku, vatn), order.by = river$Date)

# optional: difference (if needed)
# r_joku <- diff(joku)
# r_vatn <- diff(vatn)
# data_ret <- zoo::zoo(
#   cbind(r_joku, r_vatn),
#   order.by = river$Date[-1]
# )

# rename
colnames(data) <- c("ljoku", "lvatn")

# convert to plain matrix (for estimation)
mx <- zoo::coredata(data)
