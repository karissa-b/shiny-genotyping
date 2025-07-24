# genotyping full experiment fish
library(tidyverse)
library(readxl)
library(magrittr)
library(RColorBrewer)
library(scales)
library(ggfortify)

theme_set(theme_linedraw())

# read in the data
files <- list.files("~/Downloads", pattern = "20250723AW degs2P6", full.names = TRUE)
files

# melt cureve deriviities

dfs <- read_tsv(files[1], skip = 8) %>%
  pivot_longer(names_to = "Reading",
               values_to = "dF",
               starts_with("Reading")
  ) %>%
  dplyr::select(Position = `Well Location`, Reading, dF)

# temperature readings
temps <- read_tsv(files[2], skip = 8) %>%
  pivot_longer(names_to = "Reading",
               values_to = "Temperature",
               starts_with("Reading")
  ) %>%
  dplyr::select(Position = `Well Location`, Reading, Temperature)

# join these together
data.meltcurve <-
  dfs %>%
  left_join(temps)


# plot out the genotyping results
data.meltcurve %>%
  ggplot(aes(x = Temperature, y = dF)) +
  geom_line(aes(group = Position)) +
  scale_x_continuous(limits = c(78,90),
                     breaks = seq(70,90))
  scale_y_continuous(labels = comma, limits = c(0,80000))

# define the hets and homs
## define mina dn maxs

mut_temp_min <- 79
mut_temp_max <- 81.8

wt_temp_max <- 85
wt_temp_min <-82

wt_deriv_max <- 200000
wt_deriv_min <- 110000

mut_deriv_max <- 200000
mut_deriv_min <- 110000


peak_mut_positions <-
  data.meltcurve %>%
  group_by(Position) %>%
  dplyr::filter(dF %>% between(mut_deriv_min, mut_deriv_max)) %>%
  dplyr::filter(Temperature %>% between(mut_temp_min, mut_temp_max)) %>%
  dplyr::filter(lag(dF, 1) < dF & lead(dF, 1) < dF) %>%
  select(Temperature, dF) %>%
  .$Position %>% unique

peak_wt_positions <-
  data.meltcurve %>%
  group_by(Position) %>%
  dplyr::filter(dF %>% between(wt_deriv_min, wt_deriv_max)) %>%
  dplyr::filter(Temperature %>% between(wt_temp_min , wt_temp_max)) %>%
  dplyr::filter(lag(dF, 1) < dF & lead(dF, 1) < dF) %>%
  select(Temperature, dF) %>%
  .$Position %>% unique

# make a genotype data frame
genotypes <- tibble(

  Position = data.meltcurve$Position %>% unique,

  wt_allele = case_when(
    Position %in% peak_wt_positions ~ TRUE,
    TRUE ~ FALSE),
  mut_allele = case_when(
    Position %in% peak_mut_positions ~ TRUE,
    TRUE ~ FALSE),
  genotype = case_when(
    wt_allele == FALSE & mut_allele == FALSE ~ "NA",
    wt_allele == TRUE & mut_allele == FALSE ~ "wt",
    wt_allele == TRUE & mut_allele == TRUE ~ "het",
    wt_allele == FALSE & mut_allele == TRUE ~ "hom"
  ) %>%
    factor(levels = c("wt", "het", "hom", "NA"))
)

# plot
data.meltcurve %>%
  left_join(genotypes) %>%
  # mutate(
  #   genotype = case_when(
  #     Temperature == 78 & dF > 22000 ~ "NA",
  #     TRUE ~ genotype
  #   )
  # ) %>%
  ggplot(
    aes(x = Temperature, y = dF, colour = genotype)
  ) +
  geom_line(
    aes(group =Position)
  ) +
  scale_y_continuous(labels = comma,
                     limits = c(0,200000)) +
  scale_x_continuous( # zoom into region of interest.
    limits = c(78, 87),
    breaks = 1:90
  ) +
  scale_color_manual(
    values = c("darkgreen", "red", "darkred","grey50")
  ) +
  facet_wrap(~genotype, nrow = 1) +
  # mut allele
  annotate(
    geom = "rect",
    alpha = 0.5,
    fill = "red",
    xmin = mut_temp_min,
    xmax = mut_temp_max,
    ymin = mut_deriv_min,
    ymax = mut_deriv_max
  ) +
  # wt allele
  annotate(
    geom = "rect",
    alpha = 0.5,
    fill = "green",
    xmin = wt_temp_min,
    xmax = wt_temp_max,
    ymin = wt_deriv_min,
    ymax = wt_deriv_max
  ) +
  theme(text = element_text(size = 18))

# unsure samples
#c(B10, C9, C10, D12, F11, )

data.meltcurve %>%
  left_join(genotypes) %>%
  #dplyr::filter(grepl(Position, pattern = "^B")) %>%
  ggplot(
    aes(x = Temperature, y = dF, colour = Position)
  ) +
  geom_line(
    aes(group = Position)
  ) +
  geom_vline(xintercept = 82) +
  #scale_color_brewer(palette = "Paired") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous( # zoom into region of interest.
    limits = c(78, 87),
    breaks = 1:90
  )

data.meltcurve %>%
  dplyr::select(-Temperature) %>%
  pivot_wider(
    names_from = "Reading", values_from = "dF",
    id_cols = Position
  ) %>%
  column_to_rownames("Position") %>%
  prcomp %>%
  autoplot(
    data = .$x %>%
      as.data.frame() %>%
      rownames_to_column("Position") %>%
      left_join(genotypes),
    colour = "genotype"
  )

meta %>%
  left_join(genotypes, by = c("pos" = "Position")) %>%
  write.csv("~/Downloads/AAgeno/nagluGenosKB.csv")

meta %<>%
  left_join(genotypes, by = c("pos" = "Position"))

read.csv("~/Downloads/AAgeno/nagluGenosKB.csv") %>%
  left_join(data.meltcurve, by = c("pos" = "Position")) %>%

#  dplyr::filter(genotype == "wt") %>%

  ggplot(
    aes(x = Temperature, y = dF, colour = repeat_naglu.)
  ) +
  geom_line(
    aes(group =pos)
  ) +
  scale_y_continuous(labels = comma,
                     limits = c(0,80000)) +
  scale_x_continuous( # zoom into region of interest.
    limits = c(78, 87),
    breaks = 1:90
  ) +
  scale_color_manual(
    values = c("darkgreen", "red", "darkred","grey50")
  ) +
  facet_wrap(~genotype)
