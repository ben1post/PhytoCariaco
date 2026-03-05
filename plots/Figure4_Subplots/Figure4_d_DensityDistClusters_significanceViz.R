library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

# ==============================================================================
# 1. SETUP: Define Colors and Ordering
# ==============================================================================

# Define your custom order based on the vector provided in your script
test_vars_finalnames <- c(
  "AMO", "MEI v.2", "Wind speed", "Precipitation", "Salinity",
  "SST", "21°C Isotherm", "NO$_3$", "PO$_4$", "SiO$_4$",
  "Primary Productivity", "Chlorophyll $a$", 
  "Diatoms", "Haptophytes", "Dinoflagellates", 
  "Cyanobacteria", "Nanoflagellates",
  "Genus Richness", "Shannon Index", "Pielou Index"
)

# Define the colors mapped explicitly to the cluster names used in logic
# Note: I am mapping the hex codes you provided to the specific group names
cluster_colors <- c(
  "Early C1" = "#2166AC",  # Dark Blue
  "C2"       = "#B2182B",  # Dark Red
  "Late C1"  = "#92C5DE"   # Light Blue
)

# ==============================================================================
# 2. DATA PROCESSING
# ==============================================================================

# We assume 'readable_table' is already in your environment from your previous code.
# We select only the comparison columns and the Variable column.
df_plot <- readable_table %>%
  select(Variable, `Early C1 vs C2`, `C2 vs Late C1`, `Early C1 vs Late C1`) %>%
  pivot_longer(
    cols = -Variable, 
    names_to = "Comparison", 
    values_to = "Sig_Value"
  )

# Create columns for the Left and Right dot identities based on the comparison name
df_plot <- df_plot %>%
  mutate(
    Left_Group = case_when(
      Comparison == "Early C1 vs C2"      ~ "Early C1",
      Comparison == "C2 vs Late C1"       ~ "C2",
      Comparison == "Early C1 vs Late C1" ~ "Early C1"
    ),
    Right_Group = case_when(
      Comparison == "Early C1 vs C2"      ~ "C2",
      Comparison == "C2 vs Late C1"       ~ "Late C1",
      Comparison == "Early C1 vs Late C1" ~ "Late C1"
    )
  )

# Filter for Significance Labels:
# We only want to display text if it contains an asterisk (*).
# If the value is a number (e.g., "0.38"), we set the label to empty string "".
df_plot <- df_plot %>%
  mutate(
    Label = ifelse(str_detect(Sig_Value, "\\*"), Sig_Value, "")
  )

# Set the Factor levels for Variable to ensure correct sorting on Y-axis
# We reverse the order so the first item in your list appears at the TOP of the plot
df_plot$Variable <- factor(df_plot$Variable, levels = rev(test_vars_finalnames))

# Set Factor levels for Comparison to maintain logical x-axis order
df_plot$Comparison <- factor(df_plot$Comparison, 
                             levels = c("Early C1 vs C2", "C2 vs Late C1", "Early C1 vs Late C1"))

# ==============================================================================
# 3. VISUALISATION
# ==============================================================================
# Assuming 'df_plot' and 'cluster_colors' are already created from the previous step

p_refined <- ggplot(df_plot, aes(x = Comparison, y = Variable)) +
  
  # 1. The Left Dot (nudged slightly left, smaller size)
  geom_point(aes(color = Left_Group), 
             position = position_nudge(x = -0.06), 
             size = 3) +
  
  # 2. The Right Dot (nudged slightly right, smaller size)
  geom_point(aes(color = Right_Group), 
             position = position_nudge(x = 0.06), 
             size = 3) +
  
  # 3. The Asterisks (nudged closer to the dots)
  geom_text(aes(label = Label), 
            position = position_nudge(x = 0.2), 
            hjust = 0, vjust = 0.75, 
            size = 5, fontface = "bold") +
  
  # 4. Styling
  scale_color_manual(values = cluster_colors, name = "Cluster") +
  scale_x_discrete(position = "top") + 
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    
    # Remove all grid lines
    panel.grid = element_blank(),
    
    legend.position = "bottom"
  )

# Display
print(p_refined)

# Save to PDF
ggsave("Cluster_Significance_Overview_Refined.pdf", plot = p_refined, width = 8, height = 10, device = cairo_pdf)





######################## v2 #################################


library(ggplot2)

library(tidyr)
library(dplyr)
library(stringr)

# ==============================================================================
# 1. SETUP: Colors & Orders (Same as before)
# ==============================================================================

test_vars_finalnames <- c(
  "AMO", "MEI v.2", "Wind speed", "Precipitation", "Salinity",
  "SST", "21°C Isotherm", "NO$_3$", "PO$_4$", "SiO$_4$",
  "Primary Productivity", "Chlorophyll $a$", 
  "Diatoms", "Haptophytes", "Dinoflagellates", 
  "Cyanobacteria", "Nanoflagellates",
  "Genus Richness", "Shannon Index", "Pielou Index"
)

cluster_colors <- c(
  "Early C1" = "#2166AC",  # Dark Blue
  "C2"       = "#B2182B",  # Dark Red
  "Late C1"  = "#92C5DE"   # Light Blue
)

# ==============================================================================
# 2. DATA PROCESSING
# ==============================================================================

# (Assuming 'readable_table' is in your environment)
df_plot <- readable_table %>%
  select(Variable, `Early C1 vs C2`, `C2 vs Late C1`, `Early C1 vs Late C1`) %>%
  pivot_longer(cols = -Variable, names_to = "Comparison", values_to = "Sig_Value") %>%
  mutate(
    Left_Group = case_when(
      Comparison == "Early C1 vs C2"      ~ "Early C1",
      Comparison == "C2 vs Late C1"       ~ "C2",
      Comparison == "Early C1 vs Late C1" ~ "Early C1"
    ),
    Right_Group = case_when(
      Comparison == "Early C1 vs C2"      ~ "C2",
      Comparison == "C2 vs Late C1"       ~ "Late C1",
      Comparison == "Early C1 vs Late C1" ~ "Late C1"
    ),
    Label = ifelse(str_detect(Sig_Value, "\\*"), Sig_Value, "")
  )

# ORDERING ADJUSTMENTS:
# 1. Variable: Keep the order, but reverse it so the first item is at the TOP of the facets
df_plot$Variable <- factor(df_plot$Variable, levels = test_vars_finalnames)

# 2. Comparison: Reverse order so "Early vs C2" is the top line in the stack
df_plot$Comparison <- factor(df_plot$Comparison, 
                             levels = rev(c("Early C1 vs C2", "C2 vs Late C1", "Early C1 vs Late C1")))

# ==============================================================================
# 3. VISUALISATION (Stacked Layout)
# ==============================================================================

p_stacked <- ggplot(df_plot, aes(y = Comparison, x = 1)) +
  
  # 1. The Left Dot
  geom_point(aes(color = Left_Group), 
             position = position_nudge(x = -0.15), 
             size = 3) +
  
  # 2. The Right Dot
  geom_point(aes(color = Right_Group), 
             position = position_nudge(x = 0.15), 
             size = 3) +
  
  # 3. The Asterisks (to the right of the dots)
  geom_text(aes(label = Label), 
            position = position_nudge(x = 0.4), 
            hjust = 0, vjust = 0.75, 
            size = 4.5, fontface = "bold") +
  
  # 4. Faceting: Group by Variable
  # switch = "y" moves the label to the left side
  facet_grid(Variable ~ ., switch = "y") +
  
  # 5. Styling
  scale_color_manual(values = cluster_colors, name = "Cluster") +
  theme_minimal() +
  theme(
    # Remove standard axes
    axis.title = element_blank(),
    axis.text.x = element_blank(), # Hide x-axis numbers
    axis.ticks = element_blank(),
    
    # Style the Comparison labels (Y-axis inside the facet)
    axis.text.y = element_text(size = 9, color = "grey30"),
    
    # Style the Variable labels (The strip text on the left)
    strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5, 
                                     face = "bold", size = 11, color = "black"),
    strip.background = element_blank(),
    
    # Spacing adjustments
    panel.grid = element_blank(),
    panel.spacing = unit(0.5, "lines"), # Space between variables
    legend.position = "bottom"
  )

# ==============================================================================
# 4. EXPORT
# ==============================================================================

print(p_stacked)

# IMPORTANT: Increase height significantly because we now have ~60 rows
ggsave("Cluster_Significance_Stacked.pdf", plot = p_stacked, width = 6, height = 15, device = cairo_pdf)