to_degrees <- function(x) {x %% 360}

bDLD <- data.frame(
  subunit = 1:26,
  protein = rep("bDLD", times = 26),
  angle = to_degrees(seq(0, 101.381887450623*25, by = 101.381887450623)),
  rise = seq(0, 5.34898272991671*25, by = 5.34898272991671)
)

DD <- data.frame(
  subunit = 1:26,
  protein = rep("MyD88-DD", times = 26),
  angle = to_degrees(seq(0, 100.126238716122*25, by = 100.126238716122)),
  rise = seq(0, 5.50906072153564*25, by = 5.50906072153564)
)

DHF91 <- data.frame(
  subunit = 1:18,
  protein = rep("DHF91", times = 18),
  angle = c(to_degrees(seq(0, 71.5864579304634*5, by = 71.5864579304634)), 
            to_degrees(seq(120, 120+71.5864579304634*5, by = 71.5864579304634)),
            to_degrees(seq(240, 240+71.5864579304634*5, by = 71.5864579304634))
            ),
  rise = rep(seq(0, 25.1025423732163*5, by = 25.1025423732163), times = 3)
)

helical_net <- rbind(bDLD, DD, DHF91)

ggplot(
  data = helical_net,
  aes(
    x = angle,
    y = rise,
    color = protein,
    shape = protein
  )
)+
  geom_point()+
  scale_x_continuous(
    breaks = c(0, 90, 180, 270, 360)
  )+
  scale_y_continuous(
    breaks = scales::breaks_width(25)
  )+
  theme_bw()
