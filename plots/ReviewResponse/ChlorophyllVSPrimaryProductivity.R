library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(scales)
library(lubridate)

# --- 1. DATA PREP (Original) ---
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds") 

CARIACO <- CARIACO %>%
  mutate(
    date = as.Date(paste(time_month, "-15", sep=""), format="%m-%Y-%d"),
    year = as.numeric(format(date, "%Y")),
    month_num = as.numeric(format(date, "%m")),
    month_label = factor(month.abb[month_num], levels = month.abb)
  )

# --- 2. DATA PREP (Cluster Groups) ---
# Update group names to match your requested legend labels exactly
GroupCluster1 <- CARIACO %>%
  filter((date >= as.Date("1996-01-01") & date <= as.Date("2003-12-31")) | 
           (date >= as.Date("2014-01-01") & date <= as.Date("2016-12-31"))) %>%
  mutate(group = "Cluster 1")

GroupCluster2 <- CARIACO %>%
  filter((date >= as.Date("2004-01-01") & date <= as.Date("2013-12-31")) | 
           (date >= as.Date("2017-06-01") & date <= as.Date("2017-12-31"))) %>%
  mutate(group = "Cluster 2")

ENV_DATA_groups <- rbind(GroupCluster1, GroupCluster2)

# Define the exact color palette requested
cluster_cols <- c("Cluster 1" = "blue", 
                  "Cluster 2" = "red")

# --- 3. CREATE TOP ROW PLOTS ---

# A1. Seasonal Boxplot (Chl)
p_chl_box <- ggplot(CARIACO, aes(x = month_label, y = Chlorophyll)) +
  geom_boxplot(fill = "#2ca25f", color = "black", alpha = 0.6, outlier.size = 0.5) +
  theme_cowplot() +
  labs(y = expression(Chl~a~(mg~m^{-3})), title = "A. Chlorophyll") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(b = 0)) + scale_y_log10()

# A2. Heatmap (Chl)
p_chl_map <- ggplot(CARIACO, aes(x = month_label, y = year, fill = Chlorophyll)) +
  geom_tile() +
  scale_fill_distiller(palette = "Greens", direction = 1, name = "Chl", trans = "log10") +
  scale_y_reverse(breaks = pretty_breaks(n=5)) +
  theme_cowplot() +
  labs(y = "Year") +
  theme(axis.title.x = element_blank(), legend.position = "bottom", legend.justification = "center",
        legend.key.width = unit(0.8, "cm"), plot.margin = margin(t = 0))

# B1. Seasonal Boxplot (PP)
p_pp_box <- ggplot(CARIACO, aes(x = month_label, y = PrimaryProductivity)) +
  geom_boxplot(fill = "#2b8cbe", color = "black", alpha = 0.6, outlier.size = 0.5) +
  theme_cowplot() +
  labs(y = expression(PP~(mg~C~m^{-2}~d^{-1})), title = "B. Primary Productivity") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(b = 0)) + scale_y_log10()

# B2. Heatmap (PP)
p_pp_map <- ggplot(CARIACO, aes(x = month_label, y = year, fill = PrimaryProductivity)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1, name = "PP", trans = "log10") +
  scale_y_reverse(breaks = pretty_breaks(n=5)) +
  theme_cowplot() +
  labs(y = NULL) +
  theme(axis.title.x = element_blank(), axis.text.y = element_blank(), legend.position = "bottom",
        legend.justification = "center", legend.key.width = unit(0.8, "cm"), plot.margin = margin(t = 0))

# --- 4. CREATE BOTTOM ROW PLOTS ---

# C. Correlation Plot
cor_test <- cor.test(log10(CARIACO$Chlorophyll), log10(CARIACO$PrimaryProductivity), use="complete.obs")
r_txt <- paste0("r = ", round(cor_test$estimate, 2))

p_corr <- ggplot(CARIACO, aes(x = Chlorophyll, y = PrimaryProductivity)) +
  geom_point(alpha = 0.3, size = 1.5, color = "grey30") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", size=0.8) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks() + 
  annotate("text", x = min(CARIACO$Chlorophyll, na.rm=T), y = max(CARIACO$PrimaryProductivity, na.rm=T), 
           label = r_txt, hjust = 0, vjust = 1, size = 4.5) +
  theme_cowplot(font_size = 12) +
  labs(x = expression(Chl~a~(mg~m^{-3})), 
       y = expression(PP~(mg~C~m^{-2}~d^{-1})),
       title = "C. Correlation")

# D. Chlorophyll Density by Cluster
chl_medians <- ENV_DATA_groups %>% group_by(group) %>% summarise(m = median(Chlorophyll, na.rm=T))

p_dens_chl <- ggplot(ENV_DATA_groups, aes(x = Chlorophyll, fill = group, color = group)) +
  geom_density(alpha = 0.4) +
  geom_rug(alpha = 0.5) +
  geom_vline(data = chl_medians, aes(xintercept = m, color = group), linetype="dashed") +
  scale_fill_manual(values = cluster_cols) +
  scale_color_manual(values = cluster_cols) +
  annotation_logticks(sides="b") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = expression(Chl~a~(mg~m^{-3})), y = "Density", title = "D. Chl a - Regime Comp") + scale_x_log10()

# E. PP Density by Cluster (WITH LEGEND)
pp_medians <- ENV_DATA_groups %>% group_by(group) %>% summarise(m = median(PrimaryProductivity, na.rm=T))

p_dens_pp <- ggplot(ENV_DATA_groups, aes(x = PrimaryProductivity, fill = group, color = group)) +
  geom_density(alpha = 0.4) +
  geom_rug(alpha = 0.5) +
  geom_vline(data = pp_medians, aes(xintercept = m, color = group), linetype="dashed") +
  scale_fill_manual(values = cluster_cols) +
  scale_color_manual(values = cluster_cols) +
  theme_cowplot(font_size = 12) +
  theme(
    legend.position = c(1, 1),        # Top-Right corner inside plot
    legend.justification = c(1, 1),   # Anchor point of legend box
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent"), 
    legend.title = element_blank(),   
    legend.text = element_text(size = 8), # Slightly smaller text to fit long labels
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(x = expression(PP~(mg~C~m^{-2}~d^{-1})), y = NULL, title = "E. PP - Regime Comp")  + scale_x_log10()

# --- 6. ALIGNMENT & ASSEMBLY ---

# Build Top Row directly without align_plots
col_chl <- plot_grid(p_chl_box, p_chl_map, ncol = 1, rel_heights = c(1, 1.5), align = "v", axis = "lr")
col_pp  <- plot_grid(p_pp_box, p_pp_map, ncol = 1, rel_heights = c(1, 1.5), align = "v", axis = "lr")
top_row <- plot_grid(col_chl, col_pp, ncol = 2)

# Build Bottom Row
bottom_row <- plot_grid(p_corr, p_dens_chl, p_dens_pp, nrow = 1, align = "h", axis = "tb")

# Final Assembly
final_fig <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1.8, 1))

# Display
ggdraw(final_fig) + theme(plot.background = element_rect(fill = "white", color = NA))
