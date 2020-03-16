library(tidyverse)
data("mtcars")
glimpse(mtcars)

# Commit - load libs/data

mtcarsGroup = mtcars %>% group_by(cyl) %>%
  mutate(mean_disp = mean(disp),
            max_mpg = max(mpg))

# Commit - add summary stats

ggplot(mtcarsGroup, aes(x = hp, y = max_mpg, color = gear, group = gear)) +
  geom_line() +
  facet_wrap(~vs)

# Commit - graph