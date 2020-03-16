library(tidyverse)
library(plotrix)
data("mtcars")
glimpse(mtcars)

mtcarsGroup = mtcars %>% group_by(cyl) %>%
  mutate(mean_disp = mean(disp),
            max_mpg = max(mpg),
         std_error = std.error(disp))

# plot max_mpg
ggplot(mtcarsGroup, aes(x = hp, y = max_mpg, color = gear, group = gear)) +
  geom_line() +
  facet_wrap(~vs)

# plog mean_disp
ggplot(mtcarsGroup) +
  geom_line(aes(x = cyl, y = mean_disp, color = carb, group = carb)) +
  geom_errorbar(aes(x = cyl, ymin = mean_disp - std_error, ymax = mean_disp + std_error, color = carb, group = carb))
  