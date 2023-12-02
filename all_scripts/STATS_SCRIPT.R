library(ggplot2)
library(ggrepel)
library(reshape2)
#Insert numbers seperated by commas into each variable column
x = c(23,30,38,40,47,50,54,55)
y = c(1.25,0.85,1.30,1.70,0.75,0.75,0.50,1)

x_bar = mean(x)
s_x = sd(x)
y_bar = mean(y)
s_y = sd(y)




#One sample T-test
t.test(x, conf.level = 0.95)
#Two sample T-test
t.test(x, y, var.equal = TRUE)
#Regression analysis
relation <- lm(y~x)

plot(y,x,col = "blue",
     abline(lm(x~y)),cex = 1.3,pch = 16,xlab = "X Variable",ylab = "Y Variable")


#ANOVA
table <- read.delim("~/Documents/R_SCRIPT_ANOVA_DOC.txt",sep = '\t')
table <- data.frame(table)

ggplot(as.data.frame(table), mapping = aes(x = choice,y=grade4)) + 
  geom_point(alpha = 3) +
  geom_hline(yintercept =0) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10)) +
  labs(x = "Distance from SUP locus", y = "log2 Fold Change")

ggplot(melt(table, id.vars = "grade4"), aes(value, variable, colour = grade4)) + 
  geom_point(size = 4)

one.way <- aov(LAB1+LAB2+LAB3+LAB4, data = table)
#model  <- lm(weight ~ group, data = table)
model <- lm(LAB1~LAB2, data=table)

#Residual plot and normal distribution plot
par(mfrow = c(2, 2))
plot(model)

residuals <- resid(model)
plot(fitted(model), residuals, abline(0,0))
qqnorm(residuals)
qqline(residuals)
summary(one.way)
