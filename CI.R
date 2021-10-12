
# library(mice)
# library(glmnet)
library(tidyverse)


my <- readxl::read_excel('D:/GT90001_长毒data.xlsx',
                         sheet = 'Sheet1'
                         )
my_cmax <- readxl::read_excel('D:/GT90001_长毒data.xlsx',
                              sheet = 'Sheet2',
                              col_names = F
                              )

my_auc <- readxl::read_excel('D:/GT90001_长毒data.xlsx', 
                              sheet = 'Sheet3',
                             col_names = F
                              )
nms <- LETTERS[1:13]

names(my_auc) <- nms
names(my_cmax) <- nms

foo <- my_auc %>% unite('key', c(B,C)) %>% select(-A) %>% 
  column_to_rownames('key') %>% t() %>% 
  as_tibble() %>% 
  mutate(day1R=Day1_GT90001C/Day1_GT90001N,
         day29R=Day29_GT90001C/Day29_GT90001N,
         day85R=Day85_GT90001C/Day85_GT90001N
         )


psych::geometric.mean(foo$day85R)
mean(foo$day85R)


m1 <- lm(foo$Day1_GT90001C~foo$Day1_GT90001N, data = foo)
summary(m1)
anova(m1)
confint(m1, level = 0.95)
residuals(m1)

shapiro.test(foo$Day85_GT90001C)
shapiro.test(foo$Day85_GT90001N)
haha <- as.data.frame(c(foo$Day85_GT90001C, foo$Day85_GT90001N))
haha$group <- c(rep('C',10), rep('N', 10))
names(haha) <- c('value','group')
haha$group <- as.factor(haha$group)

bartlett.test(value~group, data = haha)
car::leveneTest(value~group, data = haha)

# t.test默认是异方差的，并采用Welch方法矫正自由度
t.test(foo$Day85_GT90001C, foo$Day85_GT90001N,
       conf.level = 0.95,
       alternative = 'two.sided',
       paired = FALSE,
       var.equal = T
       )
t.test(value~group,
       data = haha,
       conf.level = 0.95,
       alternative = 'two.sided',
       paired = FALSE,
       var.equal = T
       )

wilcox.test(foo$Day85_GT90001C,
            foo$Day85_GT90001N,
            paired=FALSE,
            alternative="two.sided",
            conf.int = T,
            conf.level = 0.95
)

ggplot(foo, aes(x=Day1_GT90001C, y=Day1_GT90001N)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_bw()

# cmax --------------------------------------------------------------------

bar <- my_cmax %>% unite('key', c(B,C)) %>% select(-A) %>% 
  column_to_rownames('key') %>% t() %>% 
  as_tibble() %>% 
  mutate(day1R=Day1_GT90001C/Day1_GT90001N,
         day29R=Day29_GT90001C/Day29_GT90001N,
         day85R=Day85_GT90001C/Day85_GT90001N
  )

m2 <- lm(bar$Day85_GT90001C~bar$Day85_GT90001N, data = bar)
summary(m2)
anova(m2)
confint(m2, level = 0.9)
residuals(m2)

# 数据分布
shapiro.test(bar$Day85_GT90001C)
shapiro.test(bar$Day85_GT90001N)
haha <- as.data.frame(c(bar$Day85_GT90001C, bar$Day85_GT90001N))
haha$group <- c(rep('C',10), rep('N', 10))
names(haha) <- c('value','group')
haha$group <- as.factor(haha$group)

bartlett.test(value~group, data = haha)
car::leveneTest(value~group, data = haha)

t.test(bar$Day85_GT90001C, bar$Day85_GT90001N,
       conf.level = 0.95,
       alternative = 'two.sided',
       paired = FALSE,
       var.equal = F
)

wilcox.test(bar$Day85_GT90001C,
            bar$Day85_GT90001N,
            paired=FALSE,
            alternative="two.sided",
            conf.int = T,
            conf.level = 0.95,
            exact = T
            )
correlation::cor_test(data = bar, x = 'Day85_GT90001C',
                      y = 'Day85_GT90001N',
                      method = 'pearson',
                      ci=0.95
                      )

ggplot(bar, aes(x=Day1_GT90001C, y=Day1_GT90001N)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_bw()

# 中位值差值置信区间 ---------------------------------------------------------------

# Hodges–Lehmann或bootstrap估计法可以计算中位值置信区间
library(pairwiseCI)

pairwiseCI(days~1, data=metadata, alternative = "two.sided",
           conf.level = 0.95, 
           method ="Median.diff" #HL.diff
           )
Rmisc::CI(tfre$DMSO_1) # 均值
asbio::ci.median(tfre$DMSO_1, conf = 0.95)

library(rcompanion)
rcompanion::groupwiseMedian(time ~ 1,
                data       = d,
                conf       = 0.95,
                R          = 5000,
                percentile = TRUE,
                bca        = FALSE,
                basic      = FALSE,
                normal     = F,
                wilcox     = F,
                digits     = 3)

wilcox.test(days ~ 1,
            data = metadata,
            alternative="two.sided",
            correct=TRUE,
            conf.int=TRUE,
            conf.level=0.95)

library(boot)
Mboot = boot(t$time,
             function(x,i) median(x[i]),
             R=10000)

boot.ci(Mboot,
        conf = 0.95,
        type = c("norm", "basic" ,"perc", "bca")
)
# 
DescTools::MedianCI(
  Data$Likert,
  conf.level = 0.95,
  na.rm = FALSE,
  method = "exact",
  R = 10000
)

# try do it 
bootmed = apply(matrix(sample(x, rep=TRUE, 10^4*length(x)), nrow=10^4), 1, median)



# 如果数据不符合正态分布可采用BOX-COX变换



# 比值 置信区间 -----------------------------------------------------------------

# for RR, OR
library(epitools)

RRtable <- matrix(c(1017,2260,165,992),nrow = 2, ncol = 2)
# RR
epitools::riskratio.wald(RRtable)
# OR
oddsratio.wald(RRtable)


# 对于比值
library(PropCIs)
library(DescTools)

# binomial proportion
# 21 is total num
stats::binom.test(9, 31,
           0.5,
           alternative="two.sided",
           conf.level=0.95) #统计推断总是通过样本推测总体的

observed = c(9, 22)
total = sum(observed)
DescTools::BinomCI(observed, total,
        conf.level = 0.95,
        method = "clopper-pearson") # for both success and failure

PropCIs::exactci(9, 31,
        conf.level=0.95)

# multinomial proportion
observed = c(10, 9, 1, 1)

MultinomCI(observed,
           conf.level=0.95,
           method="sisonglaz")


# confidence intervals for a difference in proportions
# 21, 17分别为两组的total
PropCIs::diffscoreci(81, 100, 6, 28,
            conf.level=0.95
            )

wald2ci(81, 100, 6, 28,
        conf.level=0.95,
        adjust = "Wald"
        )






# 描述性统计 ------------------------------------------------------------------
Input = ("
Instructor  Location  Attendees
Ren         North      7
Ren         North     22
Ren         North      6
Ren         North     15
Ren         South     12
Ren         South     13
Ren         South     14
Ren         South     16
Stimpy      North     18
Stimpy      North     17
Stimpy      North     15
Stimpy      North      9
Stimpy      South     15
Stimpy      South     11
Stimpy      South     19
Stimpy      South     23
")

Data = read.table(textConnection(Input),header=TRUE)

# stand error of mean

















