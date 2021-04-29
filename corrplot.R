library(corrplot)
library(RColorBrewer)

# data must be all numeric. you can convert character categories to numbers.

data.testing <- na.omit(data[,c('var1', 'var2', 'var3', 'var4')])

M <- cor(data.testing)

res1 <- cor.mtest(data.testing, conf.level = .95)
res2 <- cor.mtest(data.testing, conf.level = .99)

corrplot(M, p.mat = res1$p, sig.level = .2, type='upper',
         bg='lightgrey', col = rev(brewer.pal(11,"RdBu")))
