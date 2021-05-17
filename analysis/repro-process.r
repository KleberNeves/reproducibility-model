setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("helper.r")
source("global.r")

data_dir = "/home/kleber/Dropbox/Scientific Research/Projects/Modelo Reprodutibilidade/Results2/0421_Part2_Sweep1/"

# To process the data by parts, allowing manual interruptions
zips = paste(data_dir, list.files(data_dir, ".zip$"), sep ="")
repdata_list = vector(mode = "list", length = length(zips))
total = length(zips)

process_df = tibble(i = 1:length(zips), z = zips, processed = F)
repdata = tibble()
load("repdata.RData")

pwalk(process_df, function (i, z, processed) {
  if (!processed) {
    print(paste0(i, "/", total))
    df = get.data.from.zip.rep(z)
    process_df[i, "processed"] <<- T
    repdata_list[[i]] <<- df
  }
})

repdata = rbind(repdata, bind_rows(repdata_list))

save(repdata, process_df, file = paste0("repdata.RData"))


# To process the data at once
# repdata = get.figure.data.rep(data_dir)

# save(repdata, file = paste0("repdata.RData"))