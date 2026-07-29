# deploy.R --- push MAP-D to shinyapps.io as a NEW app: <account>/mapd_app
#
# Deploys ONLY the three files the app needs, so the bundle stays tiny (~1 MB)
# and does not upload the raw data_packet/ or the source workbook.
#
# ONE-TIME SETUP (per machine):
#   install.packages("rsconnect")
#   # Get your token+secret from shinyapps.io -> Account -> Tokens -> Show:
#   rsconnect::setAccountInfo(name   = "<your-account>",
#                             token  = "<token>",
#                             secret = "<secret>")
#
# Then run:  Rscript deploy.R
# (Rebuild the data first if the workbook changed:  Rscript build_data.R)

library(rsconnect)

APP_NAME <- "mapd_app"

if (!file.exists("data/mapd.rds"))
  stop("data/mapd.rds missing -- run: Rscript build_data.R")

ACCOUNT <- "stejatshiny"   # target account (this machine also has 'chiragjp')

message("Deploying '", APP_NAME, "' to account '", ACCOUNT, "' ...")
rsconnect::deployApp(
  appDir      = ".",
  appName     = APP_NAME,
  account     = ACCOUNT,
  appFiles    = c("app.R", "plots.R", "data/mapd.rds"),
  forceUpdate = TRUE
)
message("Done -> https://", ACCOUNT, ".shinyapps.io/", APP_NAME, "/")
