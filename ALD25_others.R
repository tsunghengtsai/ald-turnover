# Source code to reproduce all the other analyses and figures 

# Required packages -------------------------------------------------------

library(tidyverse)
library(readxl)
library(broom)
library(patchwork)
library(ggvenn)
library(ggupset)

# Setting -----------------------------------------------------------------

theme_set(
    theme_classic() + 
        theme(legend.position = "bottom") + 
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.title = element_text(size = 16), 
            axis.text.y = element_text(size = 14),
            axis.text.x = element_text(size = 16),
            legend.text = element_text(size = 14), 
            title = element_text(size = 16)
        )
)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pal_mod <- c("darkgray", "darkred")

path_data <- "data/ALD_others.xlsx"

# Utility functions -------------------------------------------------------

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

annotate_pval <- function(pval) {
    ifelse(pval<0.001, "***p<0.001", ifelse(
        pval<0.01, paste0("**p=", specify_decimal(pval, 3)), ifelse(
            pval<0.05, paste0("*p=", specify_decimal(pval, 3)), paste0("p=", specify_decimal(pval, 3))
        ))
    )
}

fit_sep <- function(df) {
    df$peptide <- factor(df$peptide)
    init <- df |> 
        group_by(peptide) |> 
        summarise(min = min(labeling), amp = max(labeling) - min(labeling))
    
    nls(
        labeling ~ label_min[peptide] + label_amp[peptide] * (1 - exp(-k[condition] * time)), 
        data = df, 
        start = list(label_min = init$min, label_amp = init$amp, k = rep(0.02, 2))
    )
}

fit_shared <- function(df) {
    df$peptide <- factor(df$peptide)
    init <- df |> 
        group_by(peptide) |> 
        summarise(min = min(labeling), amp = max(labeling) - min(labeling))
    
    nls(
        labeling ~ label_min[peptide] + label_amp[peptide] * (1 - exp(-k * time)), 
        data = df, 
        start = list(label_min = init$min, label_amp = init$amp, k = 0.02)
    )
}

plot_hlabeling <- function(nested, lvl_peptide) {
    l_pep <- vector("list", nrow(nested))
    for (i in seq_along(l_pep)) {
        ts_data <- nested$data[[i]]
        ts_data$peptide <- factor(ts_data$peptide)
        ts_data$peptide_c <- factor(paste(ts_data$peptide, ts_data$condition, sep = "_"))
        
        ts_mod <- nested$mod_sep[[i]]
        ts_pred <- ts_data |> 
            distinct(condition, peptide, peptide_c) |> 
            mutate(t = list(tibble(time = 0:21))) |> 
            unnest(t)
        ts_pred <- ts_pred |> 
            mutate(labeling = predict(ts_mod, ts_pred))
        
        l_pep[[i]] <- bind_rows(
            ts_data |> mutate(type = "data"), 
            ts_pred |> mutate(type = "pred")
        )
    }
    ts_pep <- bind_rows(l_pep) |> 
        mutate(peptide = factor(peptide, levels = lvl_peptide))
    ts_data <- ts_pep |> filter(type == "data")
    ts_pred <- ts_pep |> filter(type == "pred")
    
    ts_data |> 
        ggplot(aes(time, labeling, color = peptide)) + 
        geom_point() + 
        geom_line(data = ts_pred) +
        facet_wrap(~ condition) + 
        scale_x_continuous(breaks = c(0, 3, 7, 12, 21)) +
        scale_color_manual(values = cbPalette[-c(5, 7)]) + 
        labs(x = "Time (day)", y = "[2H] labeling", color = "") + 
        theme(
            legend.position = "right",
            axis.title = element_text(size = 18), 
            axis.text = element_text(size = 14),
            legend.text = element_text(size = 12), 
            strip.text.x = element_text(size = 16)
        )
}

gg_bar <- function(df_data, df_sum, pal, ylab, plab, tlab = NULL, asize = NULL) {
    set.seed(1)
    ggplot(df_data, aes(item, norm_intensity, color = condition)) +
        geom_col(data = df_sum, position = position_dodge(0.8), width = 0.7, fill = "white") +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
        geom_errorbar(
            aes(ymin = norm_intensity-sd, ymax = norm_intensity+sd), data = df_sum, 
            width = 0.2, position = position_dodge(0.8)
        ) + 
        annotate(
            "text", label = plab, x = 1:length(plab), y = max(df_data$norm_intensity) * 1.1, 
            size = ifelse(is.null(asize), ifelse(length(plab) < 3, 6, 4), asize)
        ) + 
        scale_color_manual(values = pal) +
        labs(x = "", y = ylab, title = tlab, color = "")
}

# Read data ---------------------------------------------------------------

sheets <- path_data |> 
    excel_sheets() |> 
    set_names()

nonms_list <- map(sheets, read_excel, path = path_data)

lvl_condition <- c("Pair-fed", "Ethanol-fed")

# Fig. 1B -----------------------------------------------------------------

df <- nonms_list[["Fig1B"]]
names(df) <- c("condition", "time", "bwt_gain", "lwt", "r_lwt_bwt")
df$condition <- factor(df$condition, lvl_condition)

# Body weight gain
df2 <- df |> 
    select(condition, time, norm_intensity = bwt_gain)
df_sum <- df2 |> 
    group_by(condition) |> 
    summarise(sd = sd(norm_intensity), norm_intensity = mean(norm_intensity))
plab <- tidy(t.test(norm_intensity ~ condition, df2))$p.value |> annotate_pval()
gg_1b1 <- ggplot(df2, aes(condition, norm_intensity)) + 
    geom_col(data = df_sum, aes(color = condition), position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_jitter(aes(color = condition), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
    geom_errorbar(
        aes(ymin = norm_intensity-sd, ymax = norm_intensity+sd, color = condition), data = df_sum, 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate("text", label = plab, x = 1.5, y = max(df2$norm_intensity) * 1.1, size = 6) +
    scale_x_discrete(labels = NULL) + 
    scale_color_manual(values = pal_mod) + 
    guides(colour = guide_legend(nrow = 2)) + 
    labs(x = NULL, y = "Body weight gain (g)", color = "")

# Liver weight
df2 <- df |> 
    select(condition, time, norm_intensity = lwt)
df_sum <- df2 |> 
    group_by(condition) |> 
    summarise(sd = sd(norm_intensity), norm_intensity = mean(norm_intensity))
plab <- tidy(t.test(norm_intensity ~ condition, df2))$p.value |> annotate_pval()
gg_1b2 <- ggplot(df2, aes(condition, norm_intensity)) + 
    geom_col(data = df_sum, aes(color = condition), position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_jitter(aes(color = condition), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
    geom_errorbar(
        aes(ymin = norm_intensity-sd, ymax = norm_intensity+sd, color = condition), data = df_sum, 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate("text", label = plab, x = 1.5, y = max(df2$norm_intensity) * 1.1, size = 6) +
    scale_x_discrete(labels = NULL) + 
    scale_color_manual(values = pal_mod) + 
    guides(colour = guide_legend(nrow = 2)) + 
    labs(x = NULL, y = "Liver weight (g)", color = "")

# Liver-body weight ratio
df2 <- df |> 
    select(condition, time, norm_intensity = r_lwt_bwt)
df_sum <- df2 |> 
    group_by(condition) |> 
    summarise(sd = sd(norm_intensity), norm_intensity = mean(norm_intensity))
plab <- tidy(t.test(norm_intensity ~ condition, df2))$p.value |> annotate_pval()
gg_1b3 <- ggplot(df2, aes(condition, norm_intensity)) + 
    geom_col(data = df_sum, aes(color = condition), position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_jitter(aes(color = condition), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
    geom_errorbar(
        aes(ymin = norm_intensity-sd, ymax = norm_intensity+sd, color = condition), data = df_sum, 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate("text", label = plab, x = 1.5, y = max(df2$norm_intensity) * 1.1, size = 6) +
    scale_x_discrete(labels = NULL) + 
    scale_color_manual(values = pal_mod) + 
    guides(colour = guide_legend(nrow = 2)) + 
    labs(x = NULL, y = "Liver-body weight ratio", color = "")

gg_1b1 + gg_1b2 + gg_1b3 + guide_area() + 
    plot_layout(guides = 'collect')

# Fig. 1C -----------------------------------------------------------------

# Hepatic TG
df <- nonms_list[["Fig1C_1"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
df <- df |> 
    mutate(item = "Hepatic TG") |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pval <- glance(t.test(norm_intensity ~ condition, df))[["p.value"]]
gg_1c1 <- gg_bar(df, df_sum, pal_mod, "TG (mg/g liver tissue)", annotate_pval(pval))

# Lipid
df <- nonms_list[["Fig1C_2"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
df <- df |> 
    mutate(item = "Lipid Peroxidation") |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pval <- glance(t.test(norm_intensity ~ condition, df))[["p.value"]]
gg_1c2 <- gg_bar(df, df_sum, pal_mod, "MDA (μM)", annotate_pval(pval))

# ALT & AST
df <- nonms_list[["Fig1C_3"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
lvl <- c("ALT", "AST")
df <- df |> 
    mutate(item = factor(item, levels = lvl)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl, levels = lvl),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == "ALT")))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == "AST")))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))
gg_1c3 <- gg_bar(df, df_sum, pal_mod, "Units/L", pvals$plab)

# ALT/AST ratio
df <- df |> 
    pivot_wider(names_from = item, values_from = norm_intensity) |> 
    mutate(norm_intensity = ALT/AST, item = "ALT/AST") |> 
    select(-ALT, -AST)
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pval <- glance(t.test(norm_intensity ~ condition, df))[["p.value"]]
gg_1c4 <- gg_bar(df, df_sum, pal_mod, "Ratio", annotate_pval(pval))

gg_1c1 <- gg_1c1 + theme(legend.position = "none")
gg_1c2 <- gg_1c2 + theme(legend.position = "none")
gg_1c3 <- gg_1c3 + theme(legend.position = "none")
gg_1c4 <- gg_1c4 + theme(legend.position = "none")
(gg_1c1 | gg_1c2) / (gg_1c3 | gg_1c4)

# Fig. 2A -----------------------------------------------------------------

df <- bind_rows(
    nonms_list[["Fig2A_1"]] |> 
        select(protein = `Leading proteins`, peptide = Sequence, modification = Modifications) |> 
        filter(modification == "Unmodified") |> 
        mutate(group = "PF"),
    nonms_list[["Fig2A_2"]] |> 
        select(protein = `Leading proteins`, peptide = Sequence, modification = Modifications) |> 
        filter(modification == "Unmodified") |> 
        mutate(group = "EF")
)

# Peptide
d <- df |> 
    group_by(peptide) |> 
    summarise(`Pair-fed` = any(group == "PF"), `Ethanol-fed` = any(group == "EF"))
ggplot(d) + 
    geom_venn(
        aes(A = `Pair-fed`, B = `Ethanol-fed`), 
        fill_color = c("darkgray","darkred"), 
        show_percentage = FALSE, text_size = 5
    ) + 
    ggtitle("Identified peptides") + 
    theme_void() +
    theme(title = element_text(size = 20))

# Protein
d <- df |> 
    group_by(protein) |> 
    summarise(`Pair-fed` = any(group == "PF"), `Ethanol-fed` = any(group == "EF"))
ggplot(d) + 
    geom_venn(
        aes(A = `Pair-fed`, B = `Ethanol-fed`), 
        fill_color = c("darkgray","darkred"), 
        show_percentage = FALSE, text_size = 5
    ) + 
    ggtitle("Identified protein") + 
    theme_void() +
    theme(title = element_text(size = 20))

# Fig. 3A -----------------------------------------------------------------

df <- nonms_list[["Fig3A"]]
names(df) <- c("condition", "time", "norm_intensity", "item")

df <- df |> 
    mutate(condition = factor(condition, levels = lvl_condition)) |> 
    mutate(time = time / 24)

ggplot(df, aes(time, norm_intensity, color = condition)) + 
    geom_point(size = 4, alpha = 0.7) + 
    geom_line(linewidth = 1) + 
    scale_color_manual(values = pal_mod) + 
    labs(x = "Time (day)", y = "Body water labeling (%)", color = "", shape = "") + 
    scale_x_continuous(breaks = c(0, 1, 3, 7, 12, 21)) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.25),
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14)
    )

# Fig. 3C -----------------------------------------------------------------

df <- nonms_list[["Fig3C"]]
names(df) <- c("condition", "time", "labeling", "item")
df$time <- df$time / 24
df$condition <- factor(df$condition, lvl_condition)

mod_sep <- fit_sep(df |> rename(peptide = item))
mod_shared <- fit_shared(df |> rename(peptide = item))

atpa_pf <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df %>% filter(condition == "Pair-fed"), 
    start = list(label_min = 0, label_amp = 0.5, k = 0.02)
)
atpa_ef <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df %>% filter(condition == "Ethanol-fed"), 
    start = list(label_min = 0, label_amp = 0.5, k = 0.02)
)

df_time <- tibble(time = (0:504) / 24)
atpa_pred <- tibble(
    time = rep(df_time$time, 2),
    condition = factor(c(rep("Pair-fed", nrow(df_time)), rep("Ethanol-fed", nrow(df_time))), levels = lvl_condition),
    labeling = c(predict(atpa_pf, df_time), predict(atpa_ef, df_time))
) 

t1 <- specify_decimal(log(2) / tidy(atpa_pf)$estimate[3], 1)
t2 <- specify_decimal(log(2) / tidy(atpa_ef)$estimate[3], 1)

atpa_pred |> 
    ggplot(aes(time, labeling, color = condition)) + 
    geom_line(aes(group = condition), linewidth = 1.5) + 
    geom_point(data = df, size = 4, alpha = 0.6) + 
    scale_color_manual(
        values = pal_mod, 
        labels = c(bquote("Pair-fed ("~t[1/2] == .(t1)~")"), bquote("Ethanol-fed ("~ t[1/2] == .(t2) ~")"))  
    ) + 
    scale_x_continuous(breaks = c(0, 1, 3, 7, 12, 21)) +
    labs(
        x = "Time (day)", y = "Total labeling", color = "", title = "ATPA | GIRPAINVGLSVSR",
        subtitle = annotate_pval(anova(mod_shared, mod_sep)$"Pr(>F)"[2])
    ) + 
    theme_classic() + 
    theme(
        legend.position = "inside", legend.position.inside = c(0.7, 0.25),
        axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        legend.text = element_text(size = 16), title = element_text(size = 18)
    )

# Fig. 4A -----------------------------------------------------------------

df <- nonms_list[["Fig4A"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
df <- df |> 
    mutate(item = "Total acetylation") |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pval <- glance(t.test(norm_intensity ~ condition, df))[["p.value"]]

gg_bar(df, df_sum, pal_mod, "Relative expression level", annotate_pval(pval)) + 
    theme(legend.position = "top")

# Fig. 4B -----------------------------------------------------------------

df_histone <- nonms_list[["Fig4B"]]
names(df_histone) <- c("condition", "time", "norm_intensity", "item")

# Rename the non-distinguishable form that will be highlighted
df_histone <- df_histone |> 
    filter(!str_detect(item, "PTM")) |> 
    mutate(item = ifelse(
        item == "H4.K5acK12ac, H4.K5acK16ac, H4.K8acK12ac,H4.K8acK16ac", 
        "H4.K5acK12ac", 
        item
    ))

change_histone <- df_histone |> 
    nest(data  = -item) |> 
    mutate(ttest = map(data, \(x) tidy(t.test(norm_intensity ~ condition, x)))) |> 
    select(-data) |> 
    unnest(ttest) |> 
    mutate(lfc = log2(estimate1 / estimate2)) |> 
    select(item, d_mean = estimate, lfc, p.value) |> 
    mutate(direction = ifelse(d_mean > 0, "up", "down"))

change_histone <- change_histone |> 
    mutate(item = str_remove(item, regex("ac$", ignore_case = T)) |> str_replace_all("ac", "-")) |> 
    separate(item, into = c("histone", "site"), sep = "\\.") |> 
    mutate(sig = ifelse(p.value<0.001, "***", ifelse(
        p.value<0.01, "**", ifelse(p.value<0.05, "*", "")
    )))

gg_h3 <- change_histone |> 
    filter(histone == "H3") |> 
    mutate(site = factor(site, levels = c("K18", "K23", "K18-K23"))) |> 
    ggplot(aes(x = site, y = lfc)) + 
    geom_bar(aes(fill = direction), stat = "identity") + 
    geom_text(aes(label = sig), vjust = 0.01) + 
    scale_fill_manual(values = c("#0072B2", "#FF6666"), limits = c("down", "up")) + 
    axis_combmatrix(sep = "-", levels = c("K18", "K23"), ylim = c(-0.15, 0.84)) + 
    labs(x = "", y = bquote(~Log[2]~ "(EF/PF)"), title = "H3 histone") + 
    theme(legend.position = "none", axis.title.x = element_blank()) + 
    theme_combmatrix(
        combmatrix.panel.point.color.fill = "#F8A31B",
        combmatrix.panel.line.size = 0,
        combmatrix.label.make_space = FALSE
    )

gg_h4 <- change_histone |> 
    filter(histone == "H4") |> 
    mutate(site = factor(site, levels = c(
        "K5", "K8", "K12", "K16", "K5-K8", "K12-K16", "K5-K12", 
        "K5-K8-K12", "K5-K8-K16", "K5-K12-K16", "K8-K12-K16", "K5-K8-K12-K16"
    ))) |> 
    ggplot(aes(x = site, y = lfc)) + 
    geom_bar(aes(fill = direction), stat = "identity") + 
    geom_text(aes(label = sig), vjust = 0.01) + 
    scale_fill_manual(values = c("#0072B2", "#FF6666"), limits = c("down", "up")) + 
    axis_combmatrix(sep = "-", levels = c("K5", "K8", "K12", "K16"), ylim = c(-0.15, 0.84)) + 
    labs(x = "", y = bquote(~Log[2]~ "(EF/PF)"), title = "H4 histone") + 
    theme(legend.position = "none", axis.title.x = element_blank()) + 
    theme_combmatrix(
        combmatrix.panel.point.color.fill = "#F8A31B",
        combmatrix.panel.line.size = 0,
        combmatrix.label.make_space = FALSE
    )

gg_h3 + gg_h4 + plot_layout(widths = c(1, 3.5))

# Fig. 4C -----------------------------------------------------------------

# Peptide
d <- nonms_list[["Fig4C_1"]] |> 
    mutate(
        `Pair-fed` = str_detect(Group, "Pair"), 
        `Ethanol-fed` = str_detect(Group, "Ethanol")
    )
gg_4c1 <- ggplot(d) + 
    geom_venn(
        aes(A = `Pair-fed`, B = `Ethanol-fed`), 
        fill_color = c("darkgray","darkred"), 
        show_percentage = FALSE, text_size = 5
    ) + 
    ggtitle("Acetylated peptides") + 
    theme_void() +
    theme(title = element_text(size = 20))

# Protein
d <- nonms_list[["Fig4C_2"]] |> 
    mutate(
        `Pair-fed` = str_detect(Group, "Pair"), 
        `Ethanol-fed` = str_detect(Group, "Ethanol")
    )
gg_4c2 <- ggplot(d) + 
    geom_venn(
        aes(A = `Pair-fed`, B = `Ethanol-fed`), 
        fill_color = c("darkgray","darkred"), 
        show_percentage = FALSE, text_size = 5
    ) + 
    ggtitle("Acetylated proteins") + 
    theme_void() +
    theme(title = element_text(size = 20))

gg_4c1 / gg_4c2

# Fig. 5B -----------------------------------------------------------------

df <- nonms_list[["Fig5B"]]
names(df) <- c("condition", "time", "norm_intensity", "item")

lvl <- c("Chemotrypsin-like", "Trypsin-like", "Caspase-like")

df <- df |> 
    mutate(item = factor(item, levels = lvl)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl, levels = lvl),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[2])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[3])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

gg_bar(df, df_sum, pal_mod, "Relative specific activity (ΔFU/min/µg protein)", 
       pvals$plab, asize = 4) + 
    coord_flip() + theme(legend.position = "top")

# Fig. 5D -----------------------------------------------------------------

df <- nonms_list[["Fig5D"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
lvl <- df |> 
    distinct(item) |> 
    pull()
df <- df |> 
    mutate(item = factor(item, levels = lvl)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl, levels = lvl),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[2])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[3])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[4])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[5])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[6])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[7])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[8])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[9])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab, asize = 4) + 
    theme(legend.position = "none", axis.text.x = element_text(size = 16, angle = 20, hjust=0.5, vjust = 0.7))

# Fig. 5E -----------------------------------------------------------------

df <- nonms_list[["Fig5E"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
lvl <- df |> 
    distinct(item) |> 
    pull()
df <- df |> 
    mutate(item = factor(item, levels = lvl)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl, levels = lvl),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[2])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab) + 
    theme(legend.position = "top")

# Fig. 6C -----------------------------------------------------------------

df <- nonms_list[["Fig6C"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
lvl <- c("Acetate", "AcCoA", "AcCar", "H3K23Ac", "Palmitate", "Cholesterol")

df <- df |> 
    mutate(item = factor(item, levels = lvl)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))

ggplot(df, aes(time, norm_intensity, color = item)) + 
    geom_point(aes(shape = item), size = 4) + 
    geom_line(linewidth = 1) + 
    scale_color_manual(values = cbPalette[-c(5, 7)]) + 
    labs(x = "Time (hour)", y = "Heavy labeled (%)", color = "", shape = "") + 
    theme(
        legend.position = "right",
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14)
    )

# Fig. 7A -----------------------------------------------------------------

df <- nonms_list[["Fig7A_1"]]
names(df) <- c("condition", "time", "norm_intensity", "item")

lvl <- c("AcCoA", "CoA", "AcCoA/CoA")

df <- df |> 
    mutate(item = ifelse(item == "AcCoA/CoA\r\n", "AcCoA/CoA", item)) |> 
    mutate(item = factor(item, levels = lvl)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    protein = factor(lvl, levels = lvl),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[2])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[3])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab, asize = 6) + 
    theme(legend.position = "top")

df <- nonms_list[["Fig7A_2"]]
names(df) <- c("condition", "time", "labeling", "item")
coa_pf <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df %>% filter(condition == "Pair-fed"), 
    start = list(label_min = 0, label_amp = 10, k = 0.02)
)
coa_ef <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df %>% filter(condition == "Ethanol-fed"), 
    start = list(label_min = 0, label_amp = 10, k = 0.02)
)
df_time <- tibble(time = 0:72)
coa_pred <- tibble(
    time = rep(df_time$time, 2),
    condition = factor(c(rep("Pair-fed", nrow(df_time)), rep("Ethanol-fed", nrow(df_time))), levels = lvl_condition),
    labeling = c(predict(coa_pf, df_time), predict(coa_ef, df_time))
) 
t1 <- specify_decimal(log(2) / tidy(coa_pf)$estimate[3], 1)
t2 <- specify_decimal(log(2) / tidy(coa_ef)$estimate[3], 1)

coa_pred |> 
    ggplot(aes(time, labeling, color = condition)) + 
    geom_line(aes(group = condition), linewidth = 1.5) + 
    geom_point(data = df, size = 4, alpha = 0.6) + 
    scale_color_manual(
        values = pal_mod, 
        labels = c(bquote("Pair-fed ("~ t[1/2] == .(t1) ~")"), bquote("Ethanol-fed ("~ t[1/2] == .(t2) ~")"))  
    ) + 
    scale_x_continuous(breaks = c(0, 8, 24, 72)) +
    labs(x = "Time (hour)", y = "M1 Ac-CoA (%)", color = "", title = "Ac-CoA turnover") + 
    theme_classic() + 
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.695, 0.225),
        axis.title = element_text(size = 18), axis.text = element_text(size = 14),
        legend.text = element_text(size = 16), title = element_text(size = 18)
    )

# Fig. 7B -----------------------------------------------------------------

df <- nonms_list[["Fig7B_1"]]
names(df) <- c("condition", "time", "norm_intensity", "item")

lvl <- df |> distinct(item) |> pull()

df <- df |> 
    mutate(item = factor(item, levels = lvl)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl, levels = lvl),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[2])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[3])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab, asize = 6) + 
    theme(legend.position = "top")

df <- nonms_list[["Fig7B_2"]]
names(df) <- c("condition", "time", "labeling", "item")
car_pf <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df %>% filter(condition == "Pair-fed"), 
    start = list(label_min = 0, label_amp = 10, k = 0.02)
)
car_ef <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df %>% filter(condition == "Ethanol-fed"), 
    start = list(label_min = 0, label_amp = 10, k = 0.02)
)
df_time <- tibble(time = 0:72)
car_pred <- tibble(
    time = rep(df_time$time, 2),
    condition = factor(c(rep("Pair-fed", nrow(df_time)), rep("Ethanol-fed", nrow(df_time))), levels = lvl_condition),
    labeling = c(predict(car_pf, df_time), predict(car_ef, df_time))
) 
t1 <- specify_decimal(log(2) / tidy(car_pf)$estimate[3], 1)
t2 <- specify_decimal(log(2) / tidy(car_ef)$estimate[3], 1)

car_pred |> 
    ggplot(aes(time, labeling, color = condition)) + 
    geom_line(aes(group = condition), linewidth = 1.5) + 
    geom_point(data = df, size = 4, alpha = 0.6) + 
    scale_color_manual(
        values = pal_mod, 
        labels = c(bquote("Pair-fed ("~ t[1/2] == .(t1) ~")"), bquote("Ethanol-fed ("~ t[1/2] == .(t2) ~")"))  
    ) + 
    scale_x_continuous(breaks = c(0, 8, 24, 72)) +
    labs(x = "Time (hour)", y = "M1 Ac-Car (%)", color = "", title = "Ac-Carnitine turnover") + 
    # theme(legend.position.inside = c(0.7, 0.25)) + 
    theme_classic() + 
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.25),
        axis.title = element_text(size = 18), axis.text = element_text(size = 14),
        legend.text = element_text(size = 16), title = element_text(size = 18)
    )

df <- nonms_list[["Fig7B_3"]]
names(df) <- c("condition", "time", "labeling", "item")
h3_pf <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df %>% filter(condition == "Pair-fed"), 
    start = list(label_min = 0, label_amp = 10, k = 0.02)
)
h3_ef <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df %>% filter(condition == "Ethanol-fed"), 
    start = list(label_min = 0, label_amp = 10, k = 0.02)
)
df_time <- tibble(time = 0:72)
h3_pred <- tibble(
    time = rep(df_time$time, 2),
    condition = factor(c(rep("Pair-fed", nrow(df_time)), rep("Ethanol-fed", nrow(df_time))), levels = lvl_condition),
    labeling = c(predict(h3_pf, df_time), predict(h3_ef, df_time))
) 
t1 <- specify_decimal(log(2) / tidy(h3_pf)$estimate[3], 1)
t2 <- specify_decimal(log(2) / tidy(h3_ef)$estimate[3], 1)

h3_pred |> 
    ggplot(aes(time, labeling, color = condition)) + 
    geom_line(aes(group = condition), linewidth = 1.5) + 
    geom_point(data = df, size = 4, alpha = 0.6) + 
    scale_color_manual(
        values = pal_mod, 
        labels = c(bquote("Pair-fed ("~ t[1/2] == .(t1) ~")"), bquote("Ethanol-fed ("~ t[1/2] == .(t2) ~")"))  
    ) + 
    scale_x_continuous(breaks = c(0, 8, 24, 72)) +
    labs(x = "Time (hour)", y = "[2H]acetyl labeling (%)", color = "", 
         title = "H3K18K23 turnover") + 
    theme_classic() + 
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.25),
        axis.title = element_text(size = 18), axis.text = element_text(size = 14),
        legend.text = element_text(size = 16), title = element_text(size = 18)
    )

# Fig. 7C -----------------------------------------------------------------

df <- nonms_list[["Fig7C"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
df$item <- ifelse(df$item == "Crat", "CrAT", df$item)

lvl <- df |> distinct(item) |> pull()

df <- df |> 
    mutate(item = factor(item, levels = lvl)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl, levels = lvl),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[2])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[3])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab, asize = 5) + 
    theme(legend.position = "right")

# Fig. 7D -----------------------------------------------------------------

df <- nonms_list[["Fig7D"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
lvl <- c("NAD", "NADH", "NAD/NADH")
df <- df |> 
    mutate(item = factor(item, levels = lvl)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl, levels = lvl),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[2])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[3])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

df1 <- df |> filter(item != "NAD/NADH")
df_sum1 <- df_sum |> filter(item != "NAD/NADH")
pvals1 <- pvals |> filter(item != "NAD/NADH")
gg_7d1 <- gg_bar(df1, df_sum1, pal_mod, "Concentration (μM)", pvals1$plab)

df2 <- df |> filter(item == "NAD/NADH")
df_sum2 <- df_sum |> filter(item == "NAD/NADH")
pvals2 <- pvals |> filter(item == "NAD/NADH")
gg_7d2 <- gg_bar(df2, df_sum2, pal_mod, "Relative NAD/NADH level", pvals2$plab)

(gg_7d1 | gg_7d2) + plot_layout(widths = c(2, 1), guides = "collect") &
    theme(legend.position = "top")

# Fig. 7E -----------------------------------------------------------------

df <- nonms_list[["Fig7E"]]
names(df) <- c("condition", "time", "norm_intensity", "item")

lvl <- df |> distinct(item) |> pull()

df <- df |> 
    mutate(item = factor(item, levels = lvl)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl, levels = lvl),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[2])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[3])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl[4])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab, asize = 6) + 
    theme(legend.position = "top")

# Fig. 7F -----------------------------------------------------------------

df <- nonms_list[["Fig7F"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
df <- df |> 
    mutate(item = "Sirt3 activity") |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pval <- glance(t.test(norm_intensity ~ condition, df))[["p.value"]]

gg_bar(df, df_sum, pal_mod, "ΔOD/mg", annotate_pval(pval)) + 
    theme(legend.position = "top")

# Fig. 8A & S12A-C --------------------------------------------------------

combined <- nonms_list[["Fig8A"]]
names(combined) <- c("condition", "time", "norm_intensity", "item")

lvl <- combined |> 
    distinct(item) |> 
    pull()

combined <- combined |>
    mutate(item = factor(item, levels = lvl)) |>
    mutate(condition = factor(condition, levels = lvl_condition)) |> 
    group_by(item) |> 
    mutate(norm_intensity = norm_intensity / mean(norm_intensity[condition == "Pair-fed"])) |> 
    ungroup()
run_id <- tibble(
    condition = rep(c("Pair-fed", "Ethanol-fed"), each = 6),
    time = rep(c(0, 1, 3, 8, 24, 72), 2),
    run = str_c(rep(c("PF", "EF"), each = 6), rep(1:6, 2), sep = "_")
)
combined <- combined |>
    inner_join(run_id)

combined_filled <- combined |> 
    select(item, run, norm_intensity) |> 
    complete(item, run) |> 
    mutate(
        inty_fac = cut(
            norm_intensity, 
            breaks = c(0, 0.25, 0.5, 1, 1.5, 2, 4, 8, max(norm_intensity, na.rm=T)),
            labels = c("0-0.25", "0.25-0.5", "0.5-1", "1-1.5", "1.5-2", "2-4", "4-8", ">8")
        )
    )

# Heatmap
ggplot(combined_filled, aes(x = item, y = run, fill = inty_fac)) + 
    geom_tile(color = "white", linewidth = 0.2) + 
    geom_text(x = 7, y = 20, label = "FA oxidation", size = 4) + 
    geom_text(x = 15, y = 20, label = "TCA cycle", size = 4) + 
    geom_text(x = 28, y = 20, label = "Amino acid metabolism", size = 4) + 
    geom_vline(xintercept = c(12.5, 17.5), linetype = "dashed") + 
    geom_hline(yintercept = 6.5) + 
    guides(fill = guide_legend(title = "Normalized\nintensity")) +
    labs(x = "", y = "") + 
    scale_y_discrete(expand = c(0, 0)) + 
    scale_x_discrete(expand = c(0, 0)) + 
    scale_fill_manual(values = c("#4040FF", "#8080FF", "#BFBFFF", "#FFCCCC", "#FF9999", "#FF6666", "#FF3333", "#FF0000"), na.value = "grey90") + 
    theme_grey(base_size = 8) + 
    theme(
        legend.text=element_text(face="bold"),
        axis.ticks=element_line(linewidth=0.4),
        axis.text.x = element_text(angle = 90),
        panel.background = element_rect("gray90"),
        plot.background = element_rect("white"),
        panel.border = element_blank()
    )

# Fig S12A: FA oxidation
lvl1 <- lvl[1:12]
lvl2 <- lvl[13:17]
lvl3 <- lvl[18:36]

df <- combined |> 
    filter(item %in% lvl1[1:5]) |> 
    mutate(item = factor(item, levels = lvl1[1:5])) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl1[1:5], levels = lvl1[1:5]),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[2])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[3])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[4])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[5])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))
gg_s12a1 <- gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab, asize = 4) + 
    theme(legend.position = "none", axis.text.x = element_text(size = 16, angle = 45, hjust=1))

df <- combined |> 
    filter(item %in% lvl1[6:12]) |> 
    mutate(item = factor(item, levels = lvl1[6:12])) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl1[6:12], levels = lvl1[6:12]),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[6])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[7])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[8])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[9])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[10])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[11])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl1[12])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))
gg_s12a2 <- gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab, asize = 4) + 
    theme(legend.position = "none", axis.text.x = element_text(size = 16, angle = 45, hjust=1))

gg_s12a1 + gg_s12a2 + 
    plot_layout(widths = c(1, 1.4)) &
    plot_annotation(title = "FA oxidation") &
    theme(plot.title = element_text(size = 22))

# Fig S12B: TCA
df <- combined |> 
    filter(item %in% lvl2) |> 
    mutate(item = factor(item, levels = lvl2)) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl2, levels = lvl2),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl2[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl2[2])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl2[3])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl2[4])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl2[5])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab, "TCA cycle", asize = 5) + 
    theme(legend.position = "none", title = element_text(size = 20))

# Fig S12C: AA
df <- combined |> 
    filter(item %in% lvl3[1:4]) |> 
    mutate(item = factor(item, levels = lvl3[1:4])) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl3[1:4], levels = lvl3[1:4]),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[1])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[2])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[3])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[4])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))
gg_s12c1 <- gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab, asize = 5) + 
    theme(axis.title = element_text(size = 20), 
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 20))

df <- combined |> 
    filter(item %in% lvl3[5:19]) |> 
    mutate(item = factor(item, levels = lvl3[5:19])) |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pvals <- tibble(
    item = factor(lvl3[5:19], levels = lvl3[5:19]),
    pval = c(
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[5])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[6])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[7])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[8])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[9])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[10])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[11])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[12])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[13])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[14])))[["p.value"]],
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[15])))[["p.value"]], 
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[16])))[["p.value"]], 
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[17])))[["p.value"]], 
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[18])))[["p.value"]], 
        glance(t.test(norm_intensity ~ condition, df |> filter(item == lvl3[19])))[["p.value"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))
gg_s12c2 <- gg_bar(df, df_sum, pal_mod, "Relative expression level", pvals$plab, asize = 5) + 
    theme(axis.title = element_text(size = 20), 
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 20))

gg_s12c1 + gg_s12c2 + 
    plot_layout(widths = c(1.2, 3), guides = 'collect') &
    plot_annotation(title = "Amino acid metabolism") &
    theme(plot.title = element_text(size = 28), legend.text = element_text(size = 20))

# Fig. 8B -----------------------------------------------------------------

df <- nonms_list[["Fig8B"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
df <- df |> 
    mutate(item = "Palmitate") |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pval <- glance(t.test(norm_intensity ~ condition, df))[["p.value"]]

gg_bar(df, df_sum, pal_mod, "Relative expression level", annotate_pval(pval)) + 
    theme(legend.position = "top")

# Fig. 8C -----------------------------------------------------------------

df <- nonms_list[["Fig8C"]]
names(df) <- c("condition", "time", "labeling", "item")

coa_pf <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df %>% filter(condition == "Pair-fed"), 
    start = list(label_min = 0, label_amp = 0.5, k = 0.02)
)
coa_ef <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df %>% filter(condition == "Ethanol-fed"), 
    start = list(label_min = 0, label_amp = 0.5, k = 0.02)
)
df_time <- tibble(time = 0:504)
coa_pred <- tibble(
    time = rep(df_time$time, 2),
    condition = factor(c(rep("Pair-fed", nrow(df_time)), rep("Ethanol-fed", nrow(df_time))), levels = lvl_condition),
    labeling = c(predict(coa_pf, df_time), predict(coa_ef, df_time))
) 
t1 <- specify_decimal(log(2) / tidy(coa_pf)$estimate[3], 1)
t2 <- specify_decimal(log(2) / tidy(coa_ef)$estimate[3], 1)

coa_pred |> 
    ggplot(aes(time, labeling, color = condition)) + 
    geom_line(aes(group = condition), linewidth = 1.5) + 
    geom_point(data = df, size = 4, alpha = 0.6) + 
    scale_color_manual(
        values = pal_mod, 
        labels = c(bquote("Pair-fed ("~ t[1/2] == .(t1) ~")"), bquote("Ethanol-fed ("~ t[1/2] == .(t2) ~")"))  
    ) + 
    scale_x_continuous(breaks = c(0, 8, 24, 72, 168, 288, 504)) +
    labs(x = "Time (hour)", y = "Fractional synthesis", color = "", title = "Palmitate turnover") + 
    theme_classic() + 
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.695, 0.225),
        axis.title = element_text(size = 18), axis.text = element_text(size = 14),
        legend.text = element_text(size = 16), title = element_text(size = 18)
    )

# Fig. 8D -----------------------------------------------------------------

df <- nonms_list[["Fig8D"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
df <- df |> 
    mutate(item = "Palmitate") |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pval <- glance(t.test(norm_intensity ~ condition, df))[["p.value"]]

gg_bar(df, df_sum, pal_mod, "Newly synthesized", annotate_pval(pval)) + 
    theme(legend.position = "top")

# Fig. S7A ----------------------------------------------------------------

df <- nonms_list[["FigS7A"]]
names(df) <- c("condition", "time", "labeling", "peptide")

# Pair-fed
df_pep <- df |> 
    filter(condition == "Pair-fed") |> 
    mutate(
        condition = ifelse(str_detect(peptide, "\\("), "Acetylated", "Native") |> factor(levels = c("Native", "Acetylated")), 
        time = time / 24,
        protein = "ATPA",
        peptide = "STVAQLVKR"
    )

atpa_a <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df_pep %>% filter(condition == "Acetylated"), 
    start = list(label_min = 0.4, label_amp = 0.6, k = 0.1)
)
atpa_n <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df_pep %>% filter(condition == "Native"), 
    start = list(label_min = 0.4, label_amp = 0.6, k = 0.1)
)

df_time <- tibble(time = (0:504) / 24)
atpa_pred <- tibble(
    time = rep(df_time$time, 2),
    condition = factor(c(rep("Acetylated", nrow(df_time)), rep("Native", nrow(df_time))), levels = c("Native", "Acetylated")),
    labeling = c(predict(atpa_a, df_time), predict(atpa_n, df_time))
) 

t1 <- specify_decimal(log(2) / tidy(atpa_n)$estimate[3], 1)
t2 <- specify_decimal(log(2) / tidy(atpa_a)$estimate[3], 1)

atpa_pred |> 
    ggplot(aes(time, labeling, color = condition)) + 
    geom_line(aes(group = condition), linewidth = 1.5) + 
    geom_point(data = df_pep, size = 4, alpha = 0.6) + 
    scale_color_manual(
        values = pal_mod, 
        labels = c(bquote("STVAQLVKR ("~t[1/2] == .(t1)~")"), bquote("STVAQLVk(261)R ("~ t[1/2] == .(t2) ~")"))  
    ) + 
    scale_x_continuous(breaks = c(0, 1, 3, 7, 12, 21)) +
    labs(
        x = "Time (day)", y = "[2H] labeling", color = "", 
        title = str_c("Pair-fed (", annotate_pval(anova(fit_shared(df_pep), fit_sep(df_pep))$"Pr(>F)"[2]), ")")
    ) + 
    theme_classic() + 
    theme(
        legend.position = "inside", legend.position.inside = c(0.7, 0.2),
        axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        legend.text = element_text(size = 15), title = element_text(size = 18)
    )

# Ethanol-fed
df_pep <- df |> 
    filter(condition == "Ethanol-fed") |> 
    mutate(
        condition = ifelse(str_detect(peptide, "\\("), "Acetylated", "Native") |> factor(levels = c("Native", "Acetylated")), 
        time = time / 24,
        protein = "ATPA",
        peptide = "STVAQLVKR"
    )

atpa_a <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df_pep %>% filter(condition == "Acetylated"), 
    start = list(label_min = 0.4, label_amp = 0.6, k = 0.1)
)
atpa_n <- nls(
    labeling ~ label_min + label_amp * (1 - exp(-k * time)), 
    data = df_pep %>% filter(condition == "Native"), 
    start = list(label_min = 0.4, label_amp = 0.6, k = 0.1)
)

df_time <- tibble(time = (0:504) / 24)
atpa_pred <- tibble(
    time = rep(df_time$time, 2),
    condition = factor(c(rep("Acetylated", nrow(df_time)), rep("Native", nrow(df_time))), levels = c("Native", "Acetylated")),
    labeling = c(predict(atpa_a, df_time), predict(atpa_n, df_time))
) 

t1 <- specify_decimal(log(2) / tidy(atpa_n)$estimate[3], 1)
t2 <- specify_decimal(log(2) / tidy(atpa_a)$estimate[3], 1)

atpa_pred |> 
    ggplot(aes(time, labeling, color = condition)) + 
    geom_line(aes(group = condition), linewidth = 1.5) + 
    geom_point(data = df_pep, size = 4, alpha = 0.6) + 
    scale_color_manual(
        values = pal_mod, 
        labels = c(bquote("STVAQLVKR ("~t[1/2] == .(t1)~")"), bquote("STVAQLVk(261)R ("~ t[1/2] == .(t2) ~")"))  
    ) + 
    scale_x_continuous(breaks = c(0, 1, 3, 7, 12, 21)) +
    labs(
        x = "Time (day)", y = "[2H] labeling", color = "", 
        title = str_c("Ethanol-fed (", annotate_pval(anova(fit_shared(df_pep), fit_sep(df_pep))$"Pr(>F)"[2]), ")")
    ) + 
    theme_classic() + 
    theme(
        legend.position = "inside", legend.position.inside = c(0.7, 0.2),
        axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        legend.text = element_text(size = 15), title = element_text(size = 18)
    )

# Fig. S8A ----------------------------------------------------------------

df_pep <- nonms_list[["FigS8A"]]
names(df_pep) <- c("condition", "time", "labeling", "peptide")

df_pep <- df_pep |> 
    mutate(
        time = time / 24,
        condition = factor(condition, levels = lvl_condition),
        site = ifelse(peptide == "KQLATKAAR", "Native", ifelse(
            peptide == "KQLATk(23)AAR", "K23", "K18K23")),
        protein = "Histone"
    )
lvl_site <- c("Native", "K23", "K18K23")
nested_pep <- df_pep |> 
    nest(data = any_of(c("time", "labeling", "condition", "peptide")))

nested_pep$mod_sep <- lapply(nested_pep$data, function(x) {
    res <- try(fit_sep(x), silent = TRUE)
    if (inherits(res, 'try-error')) return(NULL)
    res
})
nested_pep$mod_shared <- lapply(nested_pep$data, function(x) {
    res <- try(fit_shared(x), silent = TRUE)
    if (inherits(res, 'try-error')) return(NULL)
    res
})

nested_pep <- nested_pep |> 
    filter(!map_lgl(mod_shared, is.null), !map_lgl(mod_sep, is.null)) |> 
    mutate(ano = map2(mod_shared, mod_sep, anova))
nested_pep <- nested_pep |> 
    mutate(pvalue = map_dbl(ano, \(x) pull(tidy(x), p.value)[2]))

pvals <- tibble(
    site = factor(lvl_site, levels = lvl_site),
    pval = c(
        filter(nested_pep, site == lvl_site[1])[["pvalue"]],
        filter(nested_pep, site == lvl_site[2])[["pvalue"]],
        filter(nested_pep, site == lvl_site[3])[["pvalue"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

s_pep <- df_pep |> 
    distinct(site, peptide) |> 
    mutate(site = factor(site, levels = lvl_site)) |> 
    arrange(site)
plot_hlabeling(nested_pep, s_pep$peptide) + ggtitle("H3 histone")

k_pep <- nested_pep |> 
    mutate(param = map(mod_sep, tidy)) |> 
    select(protein, site, param) |> 
    unnest(param) |> 
    filter(term %in% c("k1", "k2")) |> 
    select(protein:std.error) |> 
    mutate(term = ifelse(term == "k1", lvl_condition[1], lvl_condition[2])) |> 
    rename(condition = term) |> 
    mutate(
        site = factor(site, levels = lvl_site),
        condition = factor(condition, levels = lvl_condition)
    )

ggplot(k_pep, aes(site, estimate, color = condition)) +
    geom_col(position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_errorbar(
        aes(ymin = estimate-std.error, ymax = estimate+std.error), 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate(
        "text", label = pvals$plab, x = 1:length(lvl_site), size = 6,
        y = max(k_pep$estimate + k_pep$std.error) * 1.1
    ) + 
    scale_color_manual(values = pal_mod) +
    labs(x = "", y = bquote("Turnover rate ("~day^-1 ~")"), color = "") + 
    theme(
        legend.position = "top",
        axis.title = element_text(size = 18), axis.text = element_text(size = 14),
        legend.text = element_text(size = 16), title = element_text(size = 18)
    )

# Fig. S8B ----------------------------------------------------------------

df_pep <- nonms_list[["FigS8B"]]
names(df_pep) <- c("condition", "time", "labeling", "peptide")

df_pep <- df_pep |> 
    mutate(
        time = time / 24,
        condition = factor(condition, levels = lvl_condition),
        site = ifelse(peptide == "GKGGKGLGKGGAKR", "Native", ifelse(
            peptide == "GKGGKGLGk(12)GGAKR", "K12", ifelse(
                peptide == "GKGGKGLGk(12)GGAk(16)R", "K12K16", ifelse(
                    peptide == "GKGGk(8)GLG(12)GGAk(16)R", "K8K12K16", "K5K8K12K16"
                )))),
        protein = "Histone"
    )
lvl_site <- c("Native", "K12", "K12K16", "K8K12K16", "K5K8K12K16")
nested_pep <- df_pep |> 
    nest(data = any_of(c("time", "labeling", "condition", "peptide")))

nested_pep$mod_sep <- lapply(nested_pep$data, function(x) {
    res <- try(fit_sep(x), silent = TRUE)
    if (inherits(res, 'try-error')) return(NULL)
    res
})
nested_pep$mod_shared <- lapply(nested_pep$data, function(x) {
    res <- try(fit_shared(x), silent = TRUE)
    if (inherits(res, 'try-error')) return(NULL)
    res
})

nested_pep <- nested_pep |> 
    filter(!map_lgl(mod_shared, is.null), !map_lgl(mod_sep, is.null)) |> 
    mutate(ano = map2(mod_shared, mod_sep, anova))
nested_pep <- nested_pep |> 
    mutate(pvalue = map_dbl(ano, \(x) pull(tidy(x), p.value)[2]))

pvals <- tibble(
    site = factor(lvl_site, levels = lvl_site),
    pval = c(
        filter(nested_pep, site == lvl_site[1])[["pvalue"]],
        filter(nested_pep, site == lvl_site[2])[["pvalue"]],
        filter(nested_pep, site == lvl_site[3])[["pvalue"]],
        filter(nested_pep, site == lvl_site[4])[["pvalue"]],
        filter(nested_pep, site == lvl_site[5])[["pvalue"]]
    )
) |> 
    mutate(plab = annotate_pval(pval))

s_pep <- df_pep |> 
    distinct(site, peptide) |> 
    mutate(site = factor(site, levels = lvl_site)) |> 
    arrange(site)
plot_hlabeling(nested_pep, s_pep$peptide) + ggtitle("H4 histone") + 
    theme(legend.text = element_text(size = 12))

k_pep <- nested_pep |> 
    mutate(param = map(mod_sep, tidy)) |> 
    select(protein, site, param) |> 
    unnest(param) |> 
    filter(term %in% c("k1", "k2")) |> 
    select(protein:std.error) |> 
    mutate(term = ifelse(term == "k1", lvl_condition[1], lvl_condition[2])) |> 
    rename(condition = term) |> 
    mutate(
        site = factor(site, levels = lvl_site),
        condition = factor(condition, levels = lvl_condition)
    )

ggplot(k_pep, aes(site, estimate, color = condition)) +
    geom_col(position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_errorbar(
        aes(ymin = estimate-std.error, ymax = estimate+std.error), 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate(
        "text", label = pvals$plab, x = 1:length(lvl_site), size = 5,
        y = max(k_pep$estimate + k_pep$std.error) * 1.1
    ) + 
    scale_color_manual(values = pal_mod) +
    labs(x = "", y = bquote("Turnover rate ("~day^-1 ~")"), color = "") + 
    theme(
        legend.position = "top",
        axis.title = element_text(size = 18), axis.text = element_text(size = 14),
        legend.text = element_text(size = 16), title = element_text(size = 18),
        axis.text.x = element_text(size = 16, angle = 20, hjust=1)
    )

# Fig. S10A ---------------------------------------------------------------

df_prot <- nonms_list[["FigS10A"]]
names(df_prot) <- c("condition", "time", "norm_intensity", "item")
df_prot <- df_prot |> 
    filter(item == "Total ubiquitination") |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df_prot |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pval <- glance(t.test(norm_intensity ~ condition, df_prot))[["p.value"]]

gg_bar(df_prot, df_sum, pal_mod, "Relative expression level", annotate_pval(pval)) + 
    theme(legend.position = "top")

# Fig. S10B ---------------------------------------------------------------

lvl_grp <- c("PF", "PF+Inhibitor", "EF", "EF+Inhibitor")

# Trypsin-like
df <- nonms_list[["FigS10B_1"]]
names(df) <- c("condition", "time", "inhibitor", "REU")
df <- df |> 
    mutate(group = str_c(
        ifelse(condition == "Pair-fed", "PF", "EF"),
        ifelse(inhibitor == "With", "+Inhibitor", "")
    )) |> 
    mutate(group = factor(group, levels = lvl_grp))

gg_s10b1 <- ggplot(df, aes(time, REU, color = group, shape = group)) + 
    geom_point(alpha = 0.75, size = 2) + 
    scale_color_manual(name = "", labels = lvl_grp, values = rep(pal_mod, each = 2)) + 
    scale_shape_manual(name = "", labels = lvl_grp, values = c(19, 17, 19, 17)) + 
    labs(x = "Time (min)", title = "Trypsin-like")

# Chemotrypsin-like
df <- nonms_list[["FigS10B_2"]]
names(df) <- c("condition", "time", "inhibitor", "REU")
df <- df |> 
    mutate(group = str_c(
        ifelse(condition == "Pair-fed", "PF", "EF"),
        ifelse(inhibitor == "With", "+Inhibitor", "")
    )) |> 
    mutate(group = factor(group, levels = lvl_grp))

gg_s10b2 <- ggplot(df, aes(time, REU, color = group, shape = group)) + 
    geom_point(alpha = 0.75, size = 2) + 
    scale_color_manual(name = "", labels = lvl_grp, values = rep(pal_mod, each = 2)) + 
    scale_shape_manual(name = "", labels = lvl_grp, values = c(19, 17, 19, 17)) + 
    labs(x = "Time (min)", title = "Chemotrypsin-like")

# Caspase-like
df <- nonms_list[["FigS10B_3"]]
names(df) <- c("condition", "time", "inhibitor", "REU")
df <- df |> 
    mutate(group = str_c(
        ifelse(condition == "Pair-fed", "PF", "EF"),
        ifelse(inhibitor == "With", "+Inhibitor", "")
    )) |> 
    mutate(group = factor(group, levels = lvl_grp))

gg_s10b3 <- ggplot(df, aes(time, REU, color = group, shape = group)) + 
    geom_point(alpha = 0.75, size = 2) + 
    scale_color_manual(name = "", labels = lvl_grp, values = rep(pal_mod, each = 2)) + 
    scale_shape_manual(name = "", labels = lvl_grp, values = c(19, 17, 19, 17)) + 
    labs(x = "Time (min)", title = "Caspase-like")

gg_s10b1 + gg_s10b2 + gg_s10b3 + 
    plot_layout(guides = 'collect')

# Fig. S11 ----------------------------------------------------------------

df <- nonms_list[["FigS11"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
lvl <- c("KQLATKAAR", "KQLATk(23)AAR", "k(18)QLATk(23)AAR")

df <- df |> 
    mutate(item = factor(item, levels = lvl))

ggplot(df, aes(time, norm_intensity, color = item)) + 
    geom_point(aes(shape = item), size = 4, alpha = 0.7) + 
    geom_line(linewidth = 1) + 
    scale_color_manual(values = cbPalette[-c(5, 7)]) + 
    scale_x_continuous(breaks = c(0, 1, 2, 4, 6)) +
    labs(x = "Time (hour)", y = "M3 / M0", color = "", shape = "") + 
    theme(
        legend.position = "right",
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14)
    )

# Fig. S12D ---------------------------------------------------------------

df <- nonms_list[["FigS12D"]]
names(df) <- c("condition", "time", "norm_intensity", "item")
df <- df |> 
    mutate(item = "2HB") |> 
    mutate(condition = factor(condition, levels = lvl_condition))
df_sum <- df |> 
    group_by(item, condition) |> 
    summarise(
        sd = sd(norm_intensity),
        norm_intensity = mean(norm_intensity)
    )
pval <- glance(t.test(norm_intensity ~ condition, df))[["p.value"]]

gg_bar(df, df_sum, pal_mod, "Relative expression level", annotate_pval(pval)) + 
    theme(legend.position = "top")

