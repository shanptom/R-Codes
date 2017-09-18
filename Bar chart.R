library(ggplot2)
genus<-read.csv(file = "top_order_vki14.txt", sep = "\t")

# to change the order of group

genus$site_f = factor(genus$site, levels=c('VKI.14','VKI.5','VKI.2'))  

  
#faceting (split bar chart according to site)

p<-ggplot(genus, aes(y=value, x=group, color=group, fill=group)) + 
  geom_bar( stat="identity") +    
  facet_wrap(~site_f)

#changing the xlab orientation to 90 degrees
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#changing axis labels

p<-p+ labs(y = "Relative abundance",x ="Genus")
p

#grouped barchart

ggplot(genus, aes(fill=site, y=value, x=group)) + 
  geom_bar(position="dodge", stat="identity")

#stacked barchart
ggplot(genus, aes(fill=site, y=value, x=group)) + 
  geom_bar( stat="identity")

#area chart
ggplot(genus, aes(x=site, y=value, fill=group)) + geom_area() 

