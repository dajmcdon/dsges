load(here::here("data", "SWdataUpdated.Rdata"))
dat <- data.frame(t(y))
colnames(dat) <- c('hours worked','interest rate','inflation','output',
                'consumption','investment','wages')
dat$date <- seq(to = 2018.87, by=.25, length.out = nrow(dat))

dat %>% 
  pivot_longer(-date) %>%
  ggplot(aes(date, value, color = name)) +
  geom_line(size = 1.5) +
  theme_void() +
  theme(panel.background = element_rect(fill="#002145"),
        legend.position = "none") +
  coord_cartesian(ylim = c(-4,4)) +
  scale_x_continuous(expand = expansion(0)) +
  scale_color_viridis_d(option = "C")

ggsave(here::here("docs/gfx/cover.svg"), width = 16, height = 4.5)  
