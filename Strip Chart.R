library(ggplot2)
  
alpha= read.csv('alpha_diversity1.txt',sep = '\t')

positions <- c("VKI-14", "VKI-5", "VKI-2")

basic=ggplot(alpha, aes(Sample,Alpha.Diversity.Measure)) + geom_point() +facet_wrap(~ measure, scales = "free_y")+ scale_x_discrete(limits = positions)+ labs(y = "Relative abundance",x ="Sample")

basic
# Change background
basic + theme(strip.background = element_rect(colour = "red", fill = alpha("blue",0.2) ))

# Change the text 
basic + theme(strip.text.x = element_text(colour = "red", face = "bold", size=10, angle=30))

# Change the space between parts:
basic + theme(panel.spacing = unit(4, "lines"))

# to arrange facetes in single row
facet facet_wrap(nrow = =1)

df$basic = factor(df$basic, levels = c("Good's Coverage", "Shannon (H)", "Simpson Reciprocal (1/D)", "Faith's PD"))

ggplot(df, aes(x = x, y = y, color = Case)) + 
  geom_line() + 
  facet_wrap(~facet, nr = 2)
