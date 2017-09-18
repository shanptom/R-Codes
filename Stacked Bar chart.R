library(ggplot2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(plotly)
VKI<-read.csv(file="phylumab.csv")
positions <- c("VKI-14", "VKI-5", "VKI-2") #to arrange bar in desired order place in "scale_x_discrete(limits = positions)"
cl <- colors(distinct = TRUE) #to store randaom colours
fill <- sample(cl, 38) # to make 38 randomcolours
#scale_fill_manual(values = fill) adding new pallette to chart
p4 <- ggplot(VKI, aes(x=site, y=abundance, width=0.75)) + 
  geom_bar(aes(y = abundance, x = site, fill = phylum), data = VKI, stat="identity")+
  scale_fill_manual(values = fill)+ scale_x_discrete(limits = positions)
p4
p=ggplotly(p4) # to create an interactive version of the chart.
