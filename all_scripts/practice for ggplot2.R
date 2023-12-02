library(ggplot2)
options(scipen = 999)
data("midwest", package = "ggplot2")
midwest
gege <- ggplot(midwest, aes(x = Area, y = poptotal)) +
  geom_point() + geom_smooth(method = 'lm')
#geom_smooth with the method lm creates a line of 
#best fit. there are different geom layers
plot(gege) 
#always good to save these different things as 
#variables so we can edit them later/save them

#to change the limtits of the x and y axis, we can 
#delete points outside of a certain range
gege + xlim(c(0, 0.1)) + ylim(c(0, 1000000))
#this limits x to 0-0.1 range and y to 0-10000000
# another way of doing it is to zoom in on a given 
#area of the graph
gege2 <- gege + coord_cartesian(xlim = c(0,0.1), ylim = c(0,1000000))
#these may come in handy with data that has large
#amounts of outliers

gege2 + ggtitle("Area vs Population", subtitle = "From midwest dataset") +
  xlab ("Area")+
  ylab("Population")
#these commands are used to change the lables and titles of the graphs
#the subtitle option is very useful too

#in order to change the colour, we change the options inside the geom_point()
data("midwest", package = "ggplot2")
ggplot(midwest, aes(x = area, y = poptotal)) +
  geom_point(col="steelblue", size = 3) + geom_smooth(method = 'lm', col = 'pink') +
  ggtitle("Area vs Population", subtitle = "From midwest dataset") +
  xlab ("Area")+
  ylab("Population")

#the dataset we are using has many columns, which can be viewed
#by typing the name of the dataset in the console
#if we want to put colour the data points based on their
#number in another column, we must assign this with the aes()
data("midwest", package = "ggplot2")
ggplot(midwest, aes(x = area, y = poptotal)) +
  geom_point(aes(col=state), size = 3) + geom_smooth(method = 'lm', col = 'pink') +
  ggtitle("Area vs Population", subtitle = "From midwest dataset") +
  xlab ("Area")+
  ylab("Population")
#size, shape, colour and colour can all be used to 
#differentiate different points