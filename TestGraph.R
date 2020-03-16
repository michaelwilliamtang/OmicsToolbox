library(tidyverse)
data("mtcars")
glimpse(mtcars)

mtcarsGroup = mtcars %>% group_by(cyl) %>%
  mutate(mean_disp = mean(disp),
            max_mpg = max(mpg))

ggplot(mtcarsGroup, aes(x = hp, y = max_mpg, color = gear, group = gear)) +
  geom_line() +
  facet_wrap(~vs)
