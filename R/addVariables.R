setwd("C:/Google Drive/MBG/2018.03.12 Adaptation phenology")

folder = "data/data"
dataPath = "data"
outputFolder = "results"
filenames = dir(folder)
dir.create(outputFolder)

res = read.csv(file.path(dataPath, "meteo variables.csv"))


N = length(filenames)
dfList = vector("list", N)
for(i in 1:N){
   df = read.csv(file.path(folder, filenames[i]))
   df$sitename = as.character(df$sitename)
   df = merge(df, res, all.x=TRUE, sort=FALSE)
   df = df[,c("siteID",setdiff(colnames(df),"siteID"))]
   write.csv(df, file.path(outputFolder, filenames[i]), row.names=FALSE)
}
