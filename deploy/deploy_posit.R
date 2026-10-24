# Deploy the dashboard to shinyapps.io or Posit Connect Cloud with rsconnect.
#
# 1. Install NetFunrankR from GitHub (not from a local folder), so that
#    rsconnect records where the server can install it from:
#        remotes::install_github("Akmazad/NetFunrankR")
# 2. shinyapps.io: create a free account, copy the token from
#    Account > Tokens and run once:
#        rsconnect::setAccountInfo(name = "<account>", token = "<token>",
#                                  secret = "<secret>")
#    then run this script from the repository root:
#        Rscript deploy/deploy_posit.R shinyapps <account>
# 3. Posit Connect Cloud (deploys straight from a GitHub repository): run
#        Rscript deploy/deploy_posit.R manifest
#    commit deploy/app/manifest.json, and point Connect Cloud to
#    deploy/app/app.R in your repository.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "manifest"
appDir <- "deploy/app"

if (!requireNamespace("rsconnect", quietly = TRUE))
    install.packages("rsconnect")
src <- utils::packageDescription("NetFunrankR")$RemoteType
if (is.null(src))
    warning("NetFunrankR was not installed from GitHub; install it with ",
        "remotes::install_github(\"Akmazad/NetFunrankR\") before deploying")

if (mode == "manifest") {
    rsconnect::writeManifest(appDir = appDir)
    message("Wrote ", file.path(appDir, "manifest.json"))
} else if (mode == "shinyapps") {
    account <- if (length(args) > 1) args[2] else NULL
    rsconnect::deployApp(appDir = appDir, appName = "netfunrankr",
        account = account, server = "shinyapps.io", forceUpdate = TRUE)
} else {
    stop("mode must be 'manifest' or 'shinyapps'")
}
