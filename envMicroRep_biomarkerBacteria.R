library(phyloseq)
library(dplyr)
library(magrittr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(RColorBrewer)
library(microbiome)
library(Maaslin2)

################################################################################
# Figure 2B – Pielou’s Evenness
################################################################################
even.boost
# write.csv(even.boost, "./revision1_290925/envMicroRep_github/figure2B.csv", row.names = F)

##########################
library(glmmTMB)
library(emmeans)
library(dplyr)
library(purrr)

#-----------------------------#
# Function: Fit + Extract Coefs
#-----------------------------#
fit_glmm_extract <- function(df, df_name) {
  # Ensure factors
  df <- df %>%
    mutate(
      treatment = factor(treatment),
      replicate = factor(replicate),
      plot = factor(plot)
    )
  
  # Fit Gamma GLMM
  m_fit <- try(
    glmmTMB(
      Evenness ~ treatment + (1|plot/replicate),
      family = Gamma(link = "log"),
      data = df
    ),
    silent = TRUE
  )
  
  if (inherits(m_fit, "try-error")) {
    message("Model failed for ", df_name)
    return(NULL)
  }
  
  # Extract fixed effect estimates
  coef_df <- broom.mixed::tidy(m_fit, effects = "fixed", conf.int = TRUE) %>%
    mutate(dataset = df_name)
  
  return(coef_df)
}

#-----------------------------#
# Apply to even.boost only
#-----------------------------#
dfs <- list(even.boost = even.boost)

results_list <- imap(dfs, fit_glmm_extract)

# Combine into one results table
results_all <- bind_rows(results_list, .id = "dataset_name")

# Inspect
print(results_all)

results_all$term <- results_all$term %>%
  gsub("treatmentbiogas", "Biogas digestate", .) %>%
  gsub("treatmentcow", "Cow manure", .) %>%
  gsub("treatmentpig", "Pig slurry", .)

results_all$dataset_name <- results_all$dataset_name %>%
  gsub("even.boost", "Pielou’s evenness", .)


if(!is.null(results_all) && nrow(results_all) > 0) {
  p2 <- results_all %>%
    subset(., !term %in% '(Intercept)') %>%
    ggplot(aes(y = term, x = estimate, xmin = estimate- std.error, xmax = std.error + std.error, color = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#33A02C", linewidth = 0.5) + 
    geom_pointrange() +
    facet_grid(vars(dataset_name), scales = "free", space = "free") +
    scale_color_manual(values = c("#1F78B4", "#5C4A8E", "#F06664")) + guides(color ="none") +
    geom_text(aes(label = ifelse(p.value < 0.001, "***",
                                 ifelse(p.value < 0.01, "**",
                                        ifelse(p.value < 0.05, "*", "")))),
              color = "black", size = 3, vjust = -0.5) +
    labs(title = "Pielou’s Evenness", x = "Estimate (Log scale)", y = "Fertilization regimes") +
    theme_bw() +
    theme(
      legend.position = "none",  plot.margin = margin(5, 5, 5, 5),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
      strip.text.y = element_text(angle = 0, hjust = 0, face = "bold.italic", size = 10), # Horizontal facet labels for gtdb
      strip.text.x = element_text(face = "bold", size = 10), # Facet labels for component
      axis.text.y = element_text(size = 9, face = "italic"), # Make gtdb-term labels italic
      axis.title.x = element_text(margin = margin(t = 10), size = 11),
      axis.title.y = element_text(margin = margin(r = 10), size = 11),
      panel.spacing.y = unit(0.5, "cm") # Increase vertical spacing between gtdb facets
    ) 
  
  print(p2)
}
ggsave("pilous_gamma_log.jpg", path = "./revision1_290925/envMicroRep_github/", dpi=300,
       width = 7, height =3)




################################################################################
# Figure 3A – Genus stacked bar plots
################################################################################
ps.top33.ma
# write.csv(ps.top33.ma, "./revision1_290925/envMicroRep_github/figure3A.csv", row.names = F)

colourCount <- length(unique(ps.top33.ma$Genus))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

treatment_map <- c("control" = "CS", "biogas" = "BD", "cow" = "CM", "pig" = "PS")
ps.top33.ma$treatment <- factor(dplyr::recode(ps.top33.ma$treatment, !!!treatment_map))

p.genus <- ggplot(ps.top33.ma, aes(fill = Genus, y = Abundance, x = treatment)) +
  geom_bar(position = "fill", stat = "identity", colour = NA) +
  scale_fill_manual(values = getPalette(colourCount)) +
  labs(x = "", y = "Average relative abundance") +
  theme_bw() +
  theme(
    legend.text = element_text(face = "italic", size = 12),
    text = element_text(size=12),
    axis.text = element_text(size=12),
    axis.title = element_text(size=12),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(ncol = 1))
p.genus


################################################################################
# Figure 3B – Treatment vs metabolism
################################################################################
df33 <- read.csv("F:/impala/bacteriocin/picrustRostand/pstop33/pstop33KO.csv")
df33$treatment <- factor(df33$treatment, levels = c("control", "biogas", "cow", "pig"))
metabol <- c("Carbohydrate metabolism", "Methane metabolism",
             "Nitrogen metabolism", "Sulfur metabolism")
df.meta <- subset(df33, level3 %in% metabol)
df.meta1 <- merge(df.meta, meta.psshaare1, by = "SampleID")
colnames(df.meta1)[7] <- "abundance"
# write.csv(df.meta1, "./revision1_290925/envMicroRep_github/figure3B.csv", row.names = F)

df.meta1a <- dplyr::select(df.meta1, level3,  plot, treatment, replicate, abundance) %>%
  group_by(level3,  plot, treatment, replicate) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>% data.frame()


##########################
results_list <- list()

for (lv in unique(df.meta1$level3)) {
  message("Processing: ", lv)
  
  sub <- df.meta1a %>% filter(level3 == lv)
  
  ## --- Try Gamma GLMM ---
  m_fit <- try(
    glmmTMB(
      abundance ~ treatment + (1|plot/replicate),
      family = Gamma(link = "log"),
      data = sub
    ),
    silent = TRUE
  )
  
  if (!inherits(m_fit, "try-error")) {
    coef_df <- broom.mixed::tidy(m_fit, effects = "fixed", conf.int = TRUE) %>%
      mutate(level3 = lv, model = "Gamma")
    results_list[[lv]] <- coef_df
  } else {
    message("Gamma model failed for ", lv, " — fallback to Gaussian sqrt.")
    
    ## --- Fallback: sqrt + Gaussian LMM ---
    sub <- sub %>% mutate(abundance_sqrt = sqrt(abundance))
    
    m_fit2 <- try(
      lmer(abundance_sqrt ~ treatment + (1|plot/replicate), data = sub),
      silent = TRUE
    )
    
    if (!inherits(m_fit2, "try-error")) {
      coef_df <- broom.mixed::tidy(m_fit2, effects = "fixed", conf.int = TRUE) %>%
        mutate(level3 = lv, model = "Gaussian_sqrt")
      results_list[[lv]] <- coef_df
    } else {
      # If both models fail, still keep a placeholder row
      results_list[[lv]] <- tibble(
        term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA,
        conf.low = NA, conf.high = NA,
        level3 = lv, model = "failed"
      )
    }
  }
}

results_df <- bind_rows(results_list)

print(results_df)


results_df$term <- results_df$term %>%
  gsub("treatmentbiogas", "Biogas digestate", .) %>%
  gsub("treatmentcow", "Cow manure", .) %>%
  gsub("treatmentpig", "Pig slurry", .)

if(!is.null(results_df) && nrow(results_df) > 0) {
  p2 <- results_df %>%
    subset(., !term %in% '(Intercept)') %>%
    ggplot(aes(y = term, x = estimate, xmin = estimate- std.error, xmax = std.error + std.error, color = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#33A02C", linewidth = 0.5) + 
    geom_pointrange() +
    # facet_grid(~ level3, scales ="free") +
    facet_grid(vars(level3), scales = "free", space = "free") +
    scale_color_manual(values = c("#1F78B4", "#5C4A8E", "#F06664")) + guides(color ="none") +
    geom_text(aes(label = ifelse(p.value < 0.001, "***",
                                 ifelse(p.value < 0.01, "**",
                                        ifelse(p.value < 0.05, "*", "")))),
              color = "black", size = 3, vjust = -0.5) +
    labs(title = "Fertilization regimes vs Metabolism", x = "Estimate (Log scale)", y = "Fertilization regimes") +
    theme_bw() +
    theme(
      legend.position = "none",  plot.margin = margin(5, 5, 5, 5),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
      strip.text.y = element_text(angle = 0, hjust = 0, face = "bold.italic", size = 10), # Horizontal facet labels for gtdb
      strip.text.x = element_text(face = "bold", size = 10), # Facet labels for component
      axis.text.y = element_text(size = 9, face = "italic"), # Make gtdb-term labels italic
      axis.title.x = element_text(margin = margin(t = 10), size = 11),
      axis.title.y = element_text(margin = margin(r = 10), size = 11),
      panel.spacing.y = unit(0.5, "cm") # Increase vertical spacing between gtdb facets
    ) 
  
  print(p2)
}
ggsave("FertilizationvsMetabolism_gamma_log.jpg", path = "./revision1_290925/envMicroRep_github/", dpi=300,
       width = 6, height =8)

################################################################################
# Figure 4A – Core microbiome (Heatmap)
################################################################################
pseq.rel <- microbiome::transform(psshaare, "compositional")
pseq.core <- core(pseq.rel, detection = 0.0005, prevalence = .75)
pseq.core.gen <- aggregate_taxa(pseq.core, "Genus")

to_keep <- c("Acidibacter", "Arenimonas", "Bacillus", "Bradyrhizobium", 
             "Gaiella", "Mycobacterium", "Nocardioides", "Pseudolabrys", 
             "Pseudonocardia", "Reyranella", "Rubrobacter", "Solirubrobacter"
)

pseq.core.gen.fil <- subset_taxa(pseq.core.gen, Genus %in% to_keep)

prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)

p1 <- plot_core(pseq.core.gen.fil,
                plot.type = "heatmap",
                colours = rev(brewer.pal(5, "RdBu")),
                prevalences = prevalences,
                detections = detections,
                min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))") +
  ylab("Genus") +
  theme_bw() +
  theme(
    axis.text.y = element_text(face = "italic", size=12),
    axis.text.x = element_text(size=12),
    axis.title = element_text(size=12)
  )
p1 

################################################################################
# Figure 4B – Differential abundance analysis (Maaslin2)
################################################################################
mas2.dfa <- subset(mas2.df, pval <= 0.05)
# write.csv(mas2.dfa, "./revision1_290925/envMicroRep_github/Figure4B.csv", row.names = F)
mas2.dfa$coef <- mas2.dfa$coef*-1

# Add column with scientific notation
mas2.dfa <- mas2.dfa %>%
  mutate(
    pvalue.sci = formatC(pval, format = "e", digits = 2),
    hjust_pos = ifelse(coef < 0, -0.3, 1),
    label_expr = paste0("italic('p=')~'", pvalue.sci, "'")  # plotmath expression
  )

library(dplyr)

mas2.dfa <- mas2.dfa %>%
  mutate(color_group = case_when(
    coef > 0 ~ "pos",
    value == "biogas" ~ "biogas",    # example categories for value
    value == "cow" ~ "cow",
    value == "pig" ~ "pig"
  ))



p2 <- mas2.dfa %>%
  # subset(., !term %in% '(Intercept)') %>%
  ggplot(aes(y = reorder(feature, coef), x = coef, xmin = coef- stderr, xmax = coef + stderr, color = color_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#33A02C", linewidth = 0.5) + 
  geom_pointrange() +
  # facet_grid(~ Genus, scales ="free") +
  facet_grid(vars(treatment), scales = "free", space = "free") +
  scale_color_manual(values = c(
    "pos"   = "#33A02C",
    "biogas"  = "#1F78B4",
    "cow"= "#5C4A8E",
    "pig"   = "#F06664"
  ))  + guides(color ="none") +
  geom_text(aes(label = ifelse(pval < 0.001, "***",
                               ifelse(pval < 0.01, "**",
                                      ifelse(pval < 0.05, "*", "")))),
            color = "black", size = 3, vjust = -0.5) +
  labs(title = "Fertilization regimes vs Biomarker taxa", x = "Effect size (log fold change)", y = "Biomarker taxa") +
  theme_bw() +
  theme(
    legend.position = "none",  plot.margin = margin(5, 5, 5, 5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold.italic", size = 10), # Horizontal facet labels for gtdb
    strip.text.x = element_text(face = "bold", size = 10), # Facet labels for component
    axis.text.y = element_text(size = 9, face = "italic"), # Make gtdb-term labels italic
    axis.title.x = element_text(margin = margin(t = 10), size = 11),
    axis.title.y = element_text(margin = margin(r = 10), size = 11),
    panel.spacing.y = unit(0.5, "cm") # Increase vertical spacing between gtdb facets
  ) 

print(p2)

p1 <- plot_core(pseq.core.gen.fil,
                plot.type = "heatmap",
                colours = rev(brewer.pal(5, "RdBu")),
                prevalences = prevalences,
                detections = detections,
                min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))") +
  ylab("Genus") +
  theme_bw() +
  theme(
    legend.position = "none",  plot.margin = margin(5, 5, 5, 5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold.italic", size = 10), # Horizontal facet labels for gtdb
    strip.text.x = element_text(face = "bold", size = 10), # Facet labels for component
    axis.text.y = element_text(size = 9, face = "italic"), # Make gtdb-term labels italic
    axis.title.x = element_text(margin = margin(t = 10), size = 11),
    axis.title.y = element_text(margin = margin(r = 10), size = 11),
    panel.spacing.y = unit(0.5, "cm") # Increase vertical spacing between gtdb facets
  ) 
p1 


maaslin2.p1 <- ggarrange(p1, p2, ncol = 2, nrow = 1, align = "hv")
maaslin2.p1

ggsave(plot = maaslin2.p1, "core_maaslin2_290925.jpg",
       path = "./revision1_290925/envMicroRep_github/",
       dpi=300,
       width = 10, height =6)

################################################################################
# Figure 5 – Top 10 genera boxplots (example with Nocardioides)
################################################################################
#select biomarker genera
top10 <- c('Acidibacter', 'Bacillus', 'Bradyrhizobium', 'Gaiella', 
           'Mycobacterium', 'Pseudolabrys','Nocardioides', 'Pseudonocardia',
           'Reyranella', 'Solirubrobacter')
pseq.top10.m1 <- merge(pseq.top10.m, meta.psshaare1, by = "SampleID")
colnames(pseq.top10.m1)[4] <- "abundance"
# write.csv(pseq.top10.m1, "./revision1_290925/envMicroRep_github/figure5.csv", row.names = F)

pseq.top10.m1a <- dplyr::select(pseq.top10.m1, abundance, plot, treatment, Genus, replicate) %>%
  group_by(plot, treatment, Genus, replicate) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>% data.frame()

##########################
results_list1 <- list()

for (lv in unique(pseq.top10.m1a$Genus)) {
  message("Processing: ", lv)
  
  sub <- pseq.top10.m1a %>% filter(Genus == lv)
  
  # Transform abundance to (0,1) for beta regression
  sub <- sub %>% mutate(
    abundance_beta = (abundance * (nrow(sub) - 1) + 0.5) / nrow(sub)
  )
  
  ## --- Try Beta GLMM ---
  m_fit <- try(
    glmmTMB(
      abundance_beta ~ treatment + (1|plot/replicate),
      family = beta_family(link = "logit"),
      data = sub
    ),
    silent = TRUE
  )
  
  if (!inherits(m_fit, "try-error")) {
    # Extract fixed effects with back-transformed estimates
    coef_df <- broom.mixed::tidy(m_fit, effects = "fixed", conf.int = TRUE) %>%
      mutate(
        Genus = lv,
        model = "beta_family",
        estimate_backtrans = plogis(estimate),        # correct inverse logit
        conf.low_backtrans = plogis(conf.low),
        conf.high_backtrans = plogis(conf.high)
      )
    
    results_list1[[lv]] <- coef_df
    
  } else {
    message("Beta model failed for ", lv, " — fallback to Gaussian sqrt.")
    
    ## --- Fallback: sqrt + Gaussian LMM ---
    sub <- sub %>% mutate(abundance_sqrt = sqrt(abundance))
    
    m_fit2 <- try(
      lmer(abundance_sqrt ~ treatment + (1|plot/replicate), data = sub),
      silent = TRUE
    )
    
    if (!inherits(m_fit2, "try-error")) {
      coef_df <- broom.mixed::tidy(m_fit2, effects = "fixed", conf.int = TRUE) %>%
        mutate(
          Genus = lv,
          model = "Gaussian_sqrt",
          estimate_backtrans = estimate^2,          # back-transform
          conf.low_backtrans = conf.low^2,
          conf.high_backtrans = conf.high^2
        )
      results_list1[[lv]] <- coef_df
      
    } else {
      # If both models fail, keep placeholder
      results_list1[[lv]] <- tibble(
        term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA,
        conf.low = NA, conf.high = NA,
        Genus = lv, model = "failed",
        estimate_backtrans = NA,
        conf.low_backtrans = NA,
        conf.high_backtrans = NA
      )
    }
  }
}

# Combine all results
results_df_genus <- bind_rows(results_list1)
results_df_genus

results_df_genus$term <- results_df_genus$term %>%
  gsub("treatmentbiogas", "Biogas digestate", .) %>%
  gsub("treatmentcow", "Cow manure", .) %>%
  gsub("treatmentpig", "Pig slurry", .)

if(!is.null(results_df_genus) && nrow(results_df_genus) > 0) {
  p2 <- results_df_genus %>%
    subset(., !term %in% '(Intercept)') %>%
    ggplot(aes(y = term, x = estimate, xmin = estimate- std.error, xmax = estimate + std.error, color = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#33A02C", linewidth = 0.5) + 
    geom_pointrange() +
    # facet_grid(~ Genus, scales ="free") +
    facet_grid(vars(Genus), scales = "free", space = "free") +
    scale_color_manual(values = c("#1F78B4", "#5C4A8E", "#F06664")) + guides(color ="none") +
    geom_text(aes(label = ifelse(p.value < 0.001, "***",
                                 ifelse(p.value < 0.01, "**",
                                        ifelse(p.value < 0.05, "*", "")))),
              color = "black", size = 3, vjust = -0.5) +
    labs(title = "Fertilization regimes vs Biomarker taxa", x = "Estimate (Log scale)", y = "Fertilization regimes") +
    theme_bw() +
    theme(
      legend.position = "none",  plot.margin = margin(5, 5, 5, 5),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
      strip.text.y = element_text(angle = 0, hjust = 0, face = "bold.italic", size = 10), # Horizontal facet labels for gtdb
      strip.text.x = element_text(face = "bold", size = 10), # Facet labels for component
      axis.text.y = element_text(size = 9, face = "italic"), # Make gtdb-term labels italic
      axis.title.x = element_text(margin = margin(t = 10), size = 11),
      axis.title.y = element_text(margin = margin(r = 10), size = 11),
      panel.spacing.y = unit(0.5, "cm") # Increase vertical spacing between gtdb facets
    ) 
  
  print(p2)
}
ggsave("FertilizationvsBiomarker_beta_log.jpg", path = "./revision1_290925/envMicroRep_github/", dpi=300,
       width = 7, height =12)


################################################################################
# Figure 7 – Correlation plots between biomarker bacteria and functional diversity
################################################################################

# Step 6: Taxa Abundance Correlation
# write.csv(df.cor, "./revision1_290925/envMicroRep_github/figure7.csv", row.names = F)
df.cor1 <- subset(df.cor, df.cor$Genus %in% c("Bacillus", "Bradyrhizobium", "Nocardioides", "Solirubrobacter"))

##########################

library(nlme)
library(broom.mixed)
library(dplyr)
library(ggplot2)

#------------------------------------------------#
# 1. Fit models per Genus with nlme::lme
#------------------------------------------------#
results_list <- list()

for (g in unique(df.cor1$Genus)) {
  sub <- df.cor1 %>% filter(Genus == g)
  
  m <- try(
    lme(Shannon ~ mean_abund * treatment,
        random = ~1|plot/replicate,
        data = sub,
        control = lmeControl(opt = "optim")),
    silent = TRUE
  )
  
  if (!inherits(m, "try-error")) {
    coef_df <- broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE) %>%
      mutate(Genus = g)
    results_list[[g]] <- coef_df
  }
}

results_df <- bind_rows(results_list)

# Add significance stars
results_df <- results_df %>%
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ ""
  ))

results_df <- data.frame(results_df)


#------------------------------------------------#
# 2. Predict fitted values from each model
#------------------------------------------------#
pred_list <- list()

for (g in unique(df.cor1$Genus)) {
  sub <- df.cor1 %>% filter(Genus == g)
  
  m <- lme(Shannon ~ mean_abund * treatment,
           random = ~1|plot/replicate,
           data = sub,
           control = lmeControl(opt="optim"))
  
  newdat <- expand.grid(
    mean_abund = seq(min(sub$mean_abund), max(sub$mean_abund), length.out=100),
    treatment = unique(sub$treatment),
    plot = sub$plot[1],           # dummy values for random effects
    replicate = sub$replicate[1]
  )
  
  newdat$Shannon_pred <- predict(m, newdat, level=0)
  newdat$Genus <- g
  pred_list[[g]] <- newdat
}

pred_df <- bind_rows(pred_list)

#------------------------------------------------#
# 3. Plot points + model-based regression lines
#------------------------------------------------#
ggplot(df.cor1, aes(mean_abund, Shannon)) +
  geom_point(aes(colour = treatment), alpha = 0.6, size = 3) +
  geom_line(data = pred_df, aes(x = mean_abund, y = Shannon_pred, colour = treatment), linewidth=1) +
  facet_wrap(~Genus, scales = "free") +
  scale_colour_manual(values = c("#33A02C", "#1F78B4", "#5C4A8E", "#F06664"),
                      labels = c("CS", "BD", "CM", "PS")) +
  theme_bw() +
  labs(colour= "Fertilization regime", 
       x = "Mean Abundance", 
       y = "Functional (KO) diversity") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(axis.text = element_text(size=11),
        legend.position="bottom",
        plot.margin = margin(5, 5, 5, 5),
        legend.text = element_text(size=12),
        strip.text = element_text(face="italic", size=12))

##########################

cor_results <- df.cor1 %>%
  group_by(Genus) %>%
  summarise(
    cor_test = list(cor.test(Shannon, mean_abund, method = "spearman")),
    .groups = "drop"
  ) %>%
  mutate(
    estimate = map_dbl(cor_test, "estimate"),
    p.value  = map_dbl(cor_test, "p.value"),
    label = paste0("ρ = ", round(estimate, 2), 
                   ", p = ", signif(p.value, 2))
  )

cor_results$x <- cor_results$p.value
cor_results$x[1] <- 80
cor_results$x[2] <- 210
cor_results$x[3] <- 18
cor_results$x[4] <- 22

corplot <- ggplot(df.cor1, aes(mean_abund, Shannon)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(aes(colour = treatment),  size = 3) +
  
  facet_wrap(~Genus, scales = "free") +
  scale_colour_manual(values = c("#33A02C", "#1F78B4", "#5C4A8E", "#F06664"),
                      labels = c("CS", "BD", "CM", "PS")) +
  geom_text(data = cor_results, 
            aes(x = x, y = 7.81, label = label), size=4) +
  # hjust = 1.1, vjust = 1.1, inherit.aes = FALSE) +
  theme_bw() +
  labs(colour= "Fertilization regime", 
       x = "Mean Abundance", 
       y = "Functional (KO) diversity") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(axis.text = element_text(size=11),
        legend.position="bottom",
        plot.margin = margin(5, 5, 5, 5),
        legend.text = element_text(size=12),
        strip.text = element_text(face="italic", size=12))

corplot
