#store packages
tmp = installed.packages()
installedpackages = as.vector(tmp[is.na(tmp[,"Priority"]), 1])
save(installedpackages, file="~/installed_packages.rda")


#load stored package
#download 'installed_packages.rda' for fresh installations.
load("installed_packages.rda")

for (count in 1:length(installedpackages)) {
  install.packages(installedpackages[count])
}