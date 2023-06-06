library(dplyr)
library(tidyverse)

df <- read_csv("/home/gabrielsalcedo/Documentos/TYLCVD-fit/Model5/data/Resistence_interpolated_data.csv")
t <- df$Time
ip_La1582 <- df$interpolated_La1582
ip_f1tyking <- df$interpolated_f1tyking

ip_diff_La1582<-diff(ip_La1582)
ip_diff_f1tyking<-diff(ip_f1tyking)

plot(t, ip_La1582)
len =length(t)
plot(t[2:len],ip_diff_La1582, type = "o")
length(ip_La1582)
length(ip_diff_La1582)
Time <- t[2:len]
Data <- data.frame(Time = Time, incidence_La1582 = ip_diff_La1582, incidence_f1tyking = ip_diff_f1tyking)
write_csv(Data,"/home/gabrielsalcedo/Documentos/TYLCVD-fit/Model5/data/Resistence_incidence_data.csv")