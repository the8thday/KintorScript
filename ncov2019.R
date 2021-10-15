# covid-19


library(ggplot2)
library(gganimate)
library('nCov2019')


alldata <- query()

alldata$historical['China', c('jiangsu','beijing','zhejiang')] -> y
latest_data <- alldata$latest
summary(latest_data)
plot(latest_data, region = 'Global')

latest_data$table #最新文件的


p3 <- ggplot(y,
          aes(x = recovered, y=cases,  
              size = deaths,   
              colour = province
              )) + 
  geom_point(show.legend =T, alpha = 0.8) +
  scale_color_viridis_d() +
  scale_size(range = c(2, 8)) +
  #  scale_x_continuous(breaks =seq(-10,400,50))+
  # scale_x_log10() +
  labs(x = "recovered", y = "cases")


p3 + transition_time(date) +
  labs(title = "Date: {frame_time}")+
  shadow_wake(wake_length = 0.4, alpha = FALSE)




