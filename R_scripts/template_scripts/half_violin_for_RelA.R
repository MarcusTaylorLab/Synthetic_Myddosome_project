install.packages("see")
library(see)

Test_data = data.frame(
  COHORT = rep(c("3xKO", "WT", "MyD88"), each = 100),
  STIMULATION = rep(rep(c("+IL-1", "-IL-1"), each = 50), times = 3),
  Value = rep(c(rnorm(n=40, mean = 1, sd=1), rnorm(n=50, mean = .7, sd=1)), times =3)
)

color_pal <-
  c(
    "+IL-1" = "grey60", 
    "-IL-1" = "white"
  )

ggplot(
  data = Test_data,
  aes(
    x = COHORT,
    y = Value,
    fill = STIMULATION
  )
)+ 
  geom_violinhalf(
    flip = c(1,3,5),
    scale = "area",
    position = position_dodge(0.1)
  )+
  geom_boxplot(
    width = 0.2,
    outlier.shape = NA,
    alpha = 0.2,
    position = position_dodge(0.5)
  )+
  geom_pwc(
    data = Test_data,
    method = "t.test",
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")),
    label = "p.signif",
    tip.length = 0.01,
    vjust = 0.5,
    hide.ns = "p",
    dodge = 0.5
  )+
  fill_palette(
    palette = color_pal
  )
