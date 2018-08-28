# plots for montpellier

install.packages('devtools')

library(devtools)

# devtools::install_github("tidyverse/ggplot2") # install ggplot2 3.0
# this is for the labs(tag = ) functionality

library(ggplot2)
sessionInfo() # check ggplot version

library(purrr)
library(dplyr)
library(magrittr)
library(readr)

lengths <- read_csv('Desktop/Research/2018-montpellier-poster/lengths.csv') %>% 
  rename(chr = chromosome)
d <- read_csv('Desktop/Research/2018-montpellier-poster/singhaldist2k.txt')

chrom_means <- d %>% 
  select(chr, block_rate) %>%
  group_by(chr) %>% 
  summarise(mean_rho = mean(block_rate, na.rm = TRUE)) %>% 
  left_join(lengths, by = 'chr')
  
lm(mean_rho ~ lengths, data = chrom_means) %>% summary()
# R^2 = 0.48

# rho - length plot
rho_length_plot <- ggplot(chrom_means, aes(x = lengths / 1e6, y = mean_rho)) +
  geom_smooth(method = 'lm', se = FALSE, size = 2) +
  geom_point(size = 3) +
  xlab('Length (Mb)') + 
  ylab(expression(paste('Mean ', rho, 'LD (1/bp)'))) + # xlab = annotation
  theme(plot.title = element_text(family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.y = element_text(family = "Helvetica", size = 36),
        axis.title.x = element_text(family = 'Helvetica', size = 36),
        axis.text.x = element_text(family = "Helvetica", size = 36, color = 'black'),
        axis.text.y = element_text(family = "Helvetica", size = 36, color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        panel.background = element_blank(),
        plot.tag = element_text(family = "Helvetica", size = 36, color = 'black', face = 'bold')) +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  coord_cartesian(x = c(2, 10.5), y = c(0.002, 0.012)) +
  scale_y_continuous(breaks = seq(0.003, 0.012, 0.003)) +
  labs(tag = 'B') +
  annotate('text', x = 8, y = 0.011, size = 14, 
           label = 'italic(R) ^ 2 == 0.48',
           parse = TRUE) # from https://github.com/tidyverse/ggplot2/pull/1553

rho_length_plot
  
# chr 15 plot

chr15 <- filter(d, chr == 'chromosome_15')
chr15_hot <- read_csv('Desktop/Research/2018-montpellier-poster/singhaldist2k_hot_grouped.txt') %>% 
  filter(chr == 'chromosome_15') %>% 
  as_tibble() %>% 
  mutate(midpoint = (block_start + block_end) / 2) # for plot

chr15_plot <- ggplot(chr15, aes(x = (block_start / 1e6), y = block_rate)) +
  geom_line() +
  geom_point(data = chr15_hot, 
             aes(x = (midpoint / 1e6), y = block_rate + 0.06), 
             col = 'red', shape = '*', size = 14) +
  xlab('Position on chromosome 15 (Mbp)') +
  ylab(expression(paste(rho, 'LD (1/bp)'))) +
  scale_x_continuous(limits = c(0, 2.05)) +
  labs(tag = 'A') +
  theme(plot.title = element_text(family = "Helvetica", hjust = 0.5),
        axis.title.y = element_text(family = "Helvetica", size = 36, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(family = 'Helvetica', size = 36),
        axis.text.x = element_text(family = "Helvetica", size = 36, color = 'black'),
        axis.text.y = element_text(family = "Helvetica", size = 36, color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        plot.tag = element_text(family = "Helvetica", size = 36, color = 'black', face = 'bold'),
        panel.background = element_blank()) +
#  annotate('text', x = 0.3, y = 0.5, size = 6,
#           label = 'hotspots defined as \n 5-fold elevations\n relative to surrounding \n 80 kb of sequence')
  NULL

chr15_plot

# devtools::install_github("thomasp85/patchwork")
library(patchwork)

plots <- chr15_plot + rho_length_plot
ggsave('Desktop/Research/2018-montpellier-poster/combined_fig_1_labs.eps', 
       plot = plots, width = par('din')[1] * 2, height = par('din')[1])

# correlates plot

correlates <- read_csv('Desktop/Research/2018-montpellier-poster/all_correlates_cis.csv')

cols <- c('intergenic' = 'light blue', 'upstream' = 'light blue',
          'utr3' = 'dodger blue', 'CDS' = 'dodger blue', 'intronic' = 'dodger blue',
          'utr5' = 'dodger blue', 'downstream' = 'light blue', 'both' = 'light blue')

new_fixed_genplot <- ggplot(filter(correlates, correlate != 'exonic', correlate != 'is_genic'), 
                            aes(x = correlate, y = rho, fill = correlate)) + 
  geom_bar(stat = 'identity', color = 'black', size = 1) + 
  xlab('') + ylab(expression(paste('mean ', rho, 'LD (1/bp)'))) + # xlab = annotation
  geom_errorbar(aes(x = correlate, ymin = conf.low, ymax = conf.high), width = 0.2, color = 'black', size = 1) +
  scale_x_discrete(labels = c('is_genic' = 'genic', 'utr3' = "3' UTR", 'utr5' = "5' UTR"),
                   limits = c('intergenic', 'upstream', 'utr5', 'CDS', 'intronic', 'utr3', 'downstream', 'both')) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.y = element_text(family = "Helvetica", size = 32),
        axis.title.x = element_text(family = 'Helvetica', size = 32),
        axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(family = "Helvetica", size = 32, color = 'black', angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(family = "Helvetica", size = 32, color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        panel.background = element_blank(),
        axis.line.x = element_line(size = 1.2),
        axis.line.y = element_line(size = 1.2),
        axis.ticks.y = element_line(size = 1.2)) +
  geom_hline(aes(yintercept = 0.00443298965947912), linetype = 'dashed', size = 1.2) + # genomewide mean
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.0063), breaks = seq(0.000, 0.006, 0.001)) +
  guides(fill = FALSE)

new_fixed_genplot

ggsave('Desktop/Research/2018-montpellier-poster/correlates.eps', plot = new_fixed_genplot, width = par('din')[1], height = par('din')[1])

# LD rho and diversity

rho_div <- read_delim('Desktop/Research/2018-montpellier-poster/rho_diversity_genedensity_100k.txt', delim = ' ')
map_rho <- read_delim('Desktop/Research/2018-montpellier-poster/actual_map_rho_100k.txt', delim = ' ')

full <- rho_div %>%
  select(-contains('utr'), -intronic_count, -total_gene_count, -rho_count) %>%
  mutate(functional_density = CDS_count / iter_count) %>%
  left_join(map_rho %>%
              select(-contains('recombination_rho'), -iter_count, -record_count), 
              by = c('chromosome', 'start', 'end')) 


div_theme <- function(font_size = 12) {
  out <- theme(plot.title = element_text(family = "Helvetica", hjust = 0.5),
               axis.title.y = element_text(family = "Helvetica", size = font_size),
               axis.title.x = element_text(family = 'Helvetica', size = font_size),
               axis.text.x = element_text(family = "Helvetica", size = font_size, color = 'black'),
               axis.text.y = element_text(family = "Helvetica", size = font_size, color = 'black'),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               axis.line = element_line(colour = 'black', linetype = 'solid', size = 1.2),
               plot.tag = element_text(family = "Helvetica", size = font_size, color = 'black', face = 'bold'),
               panel.background = element_blank(),
               axis.ticks = element_line(size = 1.2))
  return(out)
}

rho_div_plot <- ggplot(full, aes(x = log10(rho), y = diversity)) +
  geom_point(size = 1.5) +
  div_theme(font_size = 32) +
  xlab(expression(paste('log(', rho, 'LD)'))) +
  ylab(expression(paste('Nucleotide diversity (', theta[pi], ')'))) +
  coord_cartesian(x = c(-4.5, -1)) +
  scale_x_continuous(breaks = c(-4:-1), labels = c('0.0001', 10^-3, 10^-2, 10^-1)) +
  annotate('text', x = -3.7, y = 0.043, size = 6, 
           label = 'rho == 0.609',
           parse = TRUE) +
  annotate('text', x = -3.7, y = 0.04, size = 6,
           label = 'italic(p) < 2.2 %*% 10^-16', parse = TRUE) + # plotmath
  labs(tag = 'A')
  NULL

rho_div_plot

ggsave('Desktop/Figure_3.eps', rho_div_plot, width = par('din')[1], height = par('din')[1])

ggsave('Desktop/Research/2018-montpellier-poster/rho_div.eps', plot = rho_div_plot, 
       width = par('din')[1], height = par('din')[1])

ppcor::pcor.test(full$rho, full$diversity, full$functional_density, method = 'spearman')

# R and diversity

full %<>% mutate(map_rho = 1 / (recombination_values / 1000000)) %>% # cM/Mb
  mutate(map_rho = ifelse(is.finite(map_rho), map_rho, 0))

cor.test(full$map_rho, full$diversity, method = 'spearman')

r_div_plot <- full %>% 
  mutate(map_rho = 1 / (recombination_values / 1000000)) %>% # cM/Mb
  mutate(map_rho = ifelse(is.finite(map_rho), map_rho, 0)) %>% 
  ggplot(aes(x = log10(map_rho), y = diversity)) +
  geom_point(size = 1.5) +
  div_theme(font_size = 36) +
  xlab(expression(paste('log(', R, ') (cM/Mb)'))) +
  ylab(expression(paste('Diversity (', theta[pi], ')'))) +
  coord_cartesian(x = c(-1.4, 2)) +
  scale_x_continuous(breaks = c(-1:2), labels = c(10^-1, 1, 10, 100)) +
  annotate('text', x = -0.8, y = 0.043, size = 6, 
           label = 'rho  == -0.034',
           parse = TRUE) +
  annotate('text', x = -0.8, y = 0.04, size = 6,
           label = 'italic(p) == 0.26', parse = TRUE) +
  labs(tag = 'B') 

r_div_plot

library(patchwork)
div_plots <- rho_div_plot + r_div_plot
div_plots

ggsave('Desktop/Research/2018-montpellier-poster/combined_fig_3_labs.eps', plot = div_plots, 
       width = par('din')[1] * 2, height = par('din')[1])

ppcor::pcor.test(full$map_rho, full$diversity, full$functional_density, method = 'spearman')

lm(diversity ~ map_rho, data = full) %>% summary() # nonsignificant
lm(diversity ~ rho, data = full) %>% summary()


# diversity and functional density

div_density <- full %>% 
  filter(functional_density <= 1) %>% # exclude that one bugged window
  ggplot(aes(x = functional_density, y = diversity)) +
  geom_point(size = 1.5) +
  theme(plot.title = element_text(family = "Helvetica", hjust = 0.5),
        axis.title.y = element_text(family = "Helvetica", size = 36),
        axis.title.x = element_text(family = 'Helvetica', size = 36),
        axis.text.x = element_text(family = "Helvetica", size = 36, color = 'black'),
        axis.text.y = element_text(family = "Helvetica", size = 36, color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        plot.tag = element_text(family = "Helvetica", size = 36, color = 'black', face = 'bold'),
        panel.background = element_blank()) +
  xlab('functional density (CDS / window)') +
  coord_cartesian(x = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  ylab(expression(paste('Diversity (', theta[pi], ')'))) +
  labs(tag = 'A') +
  annotate('text', x = 0.8, y = 0.013, size = 10, 
           label = 'rho  == 0.231',
           parse = TRUE) +
  annotate('text', x = 0.8, y = 0.01, size = 10,
           label = 'italic(p) < 2.2 %*% 10^-16', parse = TRUE) +
  NULL

div_density

cor.test(full$functional_density, full$diversity, method = 'spearman')

# rho and R

library(scales)

rho_r <- full %>% 
  filter(rho != 0.0, map_rho != 0.0) %>% 
  ggplot(aes(x = log10(rho), y = map_rho)) +
  geom_point(size = 1.5) +
  theme(plot.title = element_text(family = "Helvetica", hjust = 0.5),
        axis.title.y = element_text(family = "Helvetica", size = 36),
        axis.title.x = element_text(family = 'Helvetica', size = 36),
        axis.text.x = element_text(family = "Helvetica", size = 36, color = 'black'),
        axis.text.y = element_text(family = "Helvetica", size = 36, color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = 'black', linetype = 'solid'),
        plot.tag = element_text(family = "Helvetica", size = 36, color = 'black', face = 'bold'),
        panel.background = element_blank()) +
  scale_x_continuous(breaks = seq(-6, -2, 2), # ugh
                     labels = c(expression(10^-6), expression(10^-4), expression(10^-2))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(paste('log(', rho, 'LD)'))) +
  ylab(expression(paste('log(', R, ')'))) +
  labs(tag = 'B') +
  annotate('text', x = -5, y = 0.15, size = 10, 
           label = 'rho  == -0.029',
           parse = TRUE) +
  annotate('text', x = -5, y = 0.11, size = 10,
           label = 'italic(p) == 0.33', parse = TRUE) +
  NULL

rho_r

cor.test(full$rho, full$map_rho, method = 'spearman')

lm(rho ~ map_rho, data = full) %>% summary()
# R^2 = 0.003, p = 0.06

lm(diversity ~ functional_density, data = full) %>% summary()
# R^2 = 0.006, p = 0.007
cor.test(full$diversity, full$functional_density, method = 'spearman')
# spearman's rho = 0.23, p < 2.2e16

fig_4 <- div_density + rho_r
fig_4

ggsave('Desktop/Research/2018-montpellier-poster/combined_fig_4_labs.eps', plot = fig_4, 
       width = par('din')[1] * 2, height = par('din')[1])
