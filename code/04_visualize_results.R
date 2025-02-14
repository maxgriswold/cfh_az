rm(list = ls())

library(data.table)
library(ggplot2)

df <- fread("./data/processed/df_crime_analysis.csv")

# Plot average effects within non-treated sites
# and effects within treated sites. 

df_avg <- df[total_population >= 10000 & !is.na(burglary_total) & burglary_total != 0 & outlier == F, ]
df_avg[is.na(post_treat), post_treat := 0]
df_avg[, post_treat := as.factor(post_treat)]

# Remove tucson from 2021 (which didn't appear to report
# data to ucr). 
df_avg <- df_avg[!(year == 2021 & location == 'Tucson')]

df_avg[, burglary_count := (burglary_total*total_population)/1000]
df_avg[, burglary_count := sum(.SD$burglary_count), by = c("year", "post_treat")] 
df_avg[, pop_group := sum(.SD$total_population), by = c("year", "post_treat")]
df_avg[, burglary_total := 10000*(burglary_count/pop_group)]

df_avg[, assault_count := (assault_total*total_population)/1000]
df_avg[, assault_count := sum(.SD$assault_count), by = c("year", "post_treat")] 
df_avg[, pop_group := sum(.SD$total_population), by = c("year", "post_treat")]
df_avg[, assault_total := 10000*(assault_count/pop_group)]

df_avg <- unique(df_avg[, .(year, burglary_total, assault_total, post_treat)])
df_avg[, post_treat := ifelse(post_treat == 0, "Cities without a program", "Cities with Crime-Free Housing")]

df_avg <- melt(df_avg, id.vars = c("year", "post_treat"), value.name = "totals", variable.name= "crime_type")
df_avg[, crime_type := ifelse(crime_type == "burglary_total", "Burglaries", "Assaults")]

df_avg <- unique(df_avg)
  
p1 <- ggplot(df_avg, aes(x = year, y = totals, color = post_treat)) +
  geom_line(linewidth = 1) +
  facet_wrap(~crime_type, ncol = 2) +
  scale_x_continuous(limits = c(2010, 2020), breaks = seq(2010, 2020, 5)) +
  labs(title = "Assaults and burglaries do not decrease more in cities with crime-free housing \nthan cities without programs",
       x = "Year", 
       y = "Reported crimes per\n 10,000 people",
       color = "",
       caption = "Source: FBI UCR, 2010 - 2020. \nNotes: Burglaries include forced entry, non-forced entry, and attempted. \nAssaults include those with a gun, with a knife, with a weapon, \nunarmed, and simple assaults. Data displayed for cities in AZ\n w/ populations greater than 10,000.") +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(family = 'sans', size = 10),
        plot.title = element_text(family = 'sans', size = 14),
        axis.title = element_text(family = 'sans', size = 12),
        strip.text = element_text(size = 12, family = 'sans'))

ggsave(plot = p1, filename = "./figs/crime_results_simpler.pdf", height = 8.27, width = 8.27)

