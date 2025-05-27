# Source code to reproduce all the analyses and figures based on LFQ and 
# dynamic proteomics data

# Required packages -------------------------------------------------------

library(tidyverse)
library(readxl)
library(broom)
library(car)
library(patchwork)
library(ggrepel)

# Setting -----------------------------------------------------------------

theme_set(
    theme_bw() + 
        theme(
            legend.position = "bottom",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()
        )
)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pal_mod <- c("darkgray", "darkred")

# Utility functions -------------------------------------------------------

# Fit one-compartment model with separate turnover rates
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

# Fit one-compartment model with shared turnover rate
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

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall = k))

annotate_pval <- function(pval) {
    ifelse(pval<0.001, "***p<0.001", ifelse(
        pval<0.01, paste0("**p=", specify_decimal(pval, 3)), ifelse(
            pval<0.05, paste0("*p=", specify_decimal(pval, 3)), paste0("p=", specify_decimal(pval, 3))
        ))
    )
}

# LFQ protein abundance (MaxQuant) ----------------------------------------

# LFQ data of total protein abundance (Table S1A)
total <- read_excel("data/ALD_MS.xlsx", sheet = "Table S1A", range = "A2:AH1357") |> 
    select(c(
        "Majority protein IDs", "Gene names", "Protein names", paste0("PF_R", 1:9), paste0("EF_R", 1:9)
    ))
names(total) <- c("uniprot", "gene_name", "protein_name", paste0("PF_R", 1:9), paste0("EF_R", 1:9))

total <- total |> 
    mutate(across(where(is.character), ~str_replace(., "  ", " "))) |> 
    mutate(across(where(is.character), str_trim))

# nrow(total)

# Reshape from wide to long format
total_long <- total |> 
    pivot_longer(
        cols = contains("_R"), 
        names_to = "run", values_to = "intensity", 
        values_drop_na = TRUE
    ) |> 
    mutate(log2inty = log2(intensity)) |> 
    separate(run, c("condition", "rep"), remove = FALSE) |> 
    mutate(condition = factor(condition))

# Equalizing-median normalization
total_run <- total_long |>
    group_by(run) |>
    summarise(log2inty_mdn = median(log2inty, na.rm = TRUE)) |>
    mutate(log2inty_adj = median(log2inty_mdn) - log2inty_mdn)

total_long <- total_long |>
    left_join(total_run) |>
    mutate(log2inty = log2inty + log2inty_adj) |>
    select(-c(log2inty_mdn, log2inty_adj))

# ggplot(total_long, aes(run, log2inty)) + geom_boxplot(aes(fill = condition))

# Check for proteins quantified in only one run or completely missing in one group
unmatched <- total_long |> 
    group_by(uniprot) |> 
    summarise(n_cond = n_distinct(condition), .groups = "drop") |> 
    filter(n_cond != 2)

unreplicated <- total_long |> 
    group_by(uniprot, condition) |> 
    summarise(n_run = n(), .groups = "drop") |> 
    filter(n_run == 1) |> 
    distinct(uniprot)

# Nested data for significance analysis
total_nested <- total_long |> 
    anti_join(unmatched) |> 
    anti_join(unreplicated) |> 
    nest(data = any_of(c("condition", "run", "rep", "intensity", "log2inty")))

# Two-sample t-test
total_nested <- total_nested |> 
    mutate(test = map(data, ~t.test(log2inty ~ condition, data = .))) |> 
    mutate(param = map(test, tidy))

# Summarization (FC is defined as EF/PF)
total_sum <- total_nested |> 
    select(uniprot, gene_name, protein_name, param) |> 
    unnest(param) |> 
    transmute(
        uniprot = uniprot, gene_name = gene_name, protein_name = protein_name,
        log2FC = estimate, tvalue = statistic, pvalue = p.value
    ) |> 
    mutate(adj_pvalue = p.adjust(pvalue, method = "BH"))

# Significance analysis of total protein abundance (Table S1B)
total_out <- total_sum |> arrange(pvalue)

# LFQ acetylation ---------------------------------------------------------

# LFQ data of acetylated peptides (Table S4A)
aclevel <- read_excel("data/ALD_MS.xlsx", sheet = "Table S4A", range = "A2:AD174")

# Quantified values should be double
for (i in c(str_c("EF_R", 1:9), str_c("PF_R", 1:9))) {
    aclevel[[i]] <- as.double(aclevel[[i]]) 
}

# Protein information for acetylation data
annot_ac <- aclevel |>
    distinct(uniprot, swissprot, gene_name, protein_name, localization)

# Reshape from wide to long format
aclevel_long <- aclevel |> 
    select(uniprot, gene_name, peptide_idx, contains("_R")) |> 
    rename(peptide = peptide_idx) |> 
    mutate(aa_mod = as.integer(str_extract(peptide, "(?<=\\()\\d+"))) |> 
    mutate(site = str_c("K", aa_mod)) |> 
    select(-aa_mod) |> 
    pivot_longer(
        cols = contains("_R"), 
        names_to = "run", values_to = "intensity", 
        values_drop_na = TRUE
    ) |> 
    separate(run, c("condition", "rep"), remove = FALSE) |> 
    mutate(condition = factor(condition)) |> 
    filter(intensity > 0) |> 
    mutate(log2inty = log2(intensity))

# Check for acetylated peptides completely missing in one group
unmatched <- aclevel_long |> 
    group_by(uniprot, peptide) |> 
    summarise(n_cond = n_distinct(condition), .groups = "drop") |> 
    filter(n_cond != 2)

# Nested data for significance analysis
aclevel_nested <- aclevel_long |> 
    anti_join(unmatched) |> 
    nest(data = any_of(c("condition", "run", "rep", "intensity", "log2inty")))

# Two-sample t-test
aclevel_nested <- aclevel_nested |> 
    mutate(test = map(data, \(x) t.test(log2inty ~ condition, data = x))) |> 
    mutate(param = map(test, tidy))

# Summarization (FC is defined as EF/PF)
aclevel_sum <- aclevel_nested |> 
    select(uniprot, gene_name, peptide, site, param) |> 
    unnest(param) |> 
    transmute(
        uniprot = uniprot, gene_name = gene_name, site = site, peptide = peptide, 
        log2FC = estimate, tvalue = statistic, pvalue = p.value
    ) |> 
    mutate(adj_pvalue = p.adjust(pvalue, method = "BH")) |> 
    arrange(pvalue)

# Significance analysis of acetylation level (Table S4B) 
aclevel_out <- aclevel_sum |> 
    left_join(annot_ac) |>
    relocate(uniprot, swissprot, gene_name, protein_name)

# Kinetic data of native and acetylated peptides --------------------------

# Unmodified peptides (Table S2A)
unmod <- read_excel("data/ALD_MS.xlsx", sheet = "Table S2A", range = "A2:L29626")

# Acetylated peptides (Table S5A)
acetylated <- read_excel("data/ALD_MS.xlsx", sheet = "Table S5A", range = "A2:L8809")

# length(unique(unmod$uniprot))
# unmod |> distinct(uniprot, condition) |> count(uniprot) |> filter(n == 2)

# Annotate acetylated site
acetylated <- acetylated |> 
    mutate(site = map_chr(
        peptide, 
        \(x) str_c("K", unlist(str_extract_all(x, "(?<=\\()\\d+")), collapse = "")
    )) 

# Protein information
annot_turn <- bind_rows(
    acetylated |>
        distinct(uniprot, swissprot, gene_name, protein_name, localization),
    unmod |>
        distinct(uniprot, swissprot, gene_name, protein_name, localization) |>
        anti_join(acetylated |> distinct(uniprot))
)

# Nested data to assess turnover rate changes
turn_nested <- bind_rows(
    acetylated |> 
        select(uniprot, site, peptide, time, labeling, condition),
    unmod |> 
        mutate(site = "unmodified") |> 
        select(uniprot, site, peptide, time, labeling, condition)
) |> 
    mutate(condition = factor(condition)) |> 
    nest(data = any_of(c("time", "labeling", "condition", "peptide")))

# Kinetic analysis - joint modeling ---------------------------------------

# Fit models with separate and shared turnovers
turn_nested$mod_sep <- lapply(turn_nested$data, function(x) {
    res <- try(fit_sep(x), silent = TRUE)
    if (inherits(res, 'try-error')) return(NULL)
    res
})
turn_nested$mod_shared <- lapply(turn_nested$data, function(x) {
    res <- try(fit_shared(x), silent = TRUE)
    if (inherits(res, 'try-error')) return(NULL)
    res
})

# Model comparison
turn_nested <- turn_nested |> 
    filter(!map_lgl(mod_shared, is.null), !map_lgl(mod_sep, is.null)) |> 
    mutate(ano = map2(mod_shared, mod_sep, anova))

turn_nested <- turn_nested |> 
    mutate(pvalue = map_dbl(ano, \(x) pull(tidy(x), p.value)[2]))

# Inference of half-life using Delta Method
turn_nested <- turn_nested |>
    mutate(t50 = map(mod_sep, function(x) {
        bind_rows(deltaMethod(x, "log(2)/k1"), deltaMethod(x, "log(2)/k2")) |>
            select(Estimate, SE) |>
            as.data.frame() |>
            rownames_to_column("term")
    }))

# Kinetic analysis - model-based conclusions ------------------------------

# Turnover by condition
turn_k <- turn_nested |> 
    mutate(param = map(mod_sep, tidy)) |> 
    select(uniprot, site, param) |> 
    unnest(param) |> 
    filter(term %in% c("k1", "k2")) |> 
    select(uniprot:std.error) |> 
    mutate(term = ifelse(term == "k1", "EF", "PF")) |> 
    rename(condition = term)

# Half-life by condition
turn_t50 <- turn_nested |> 
    select(uniprot, site, t50) |> 
    unnest(t50) |> 
    mutate(condition = ifelse(str_detect(term, "k1"), "EF", "PF")) |> 
    select(-term)

# Turnover in wide format
turn_k_wd <- turn_k |> 
    pivot_wider(names_from = condition, values_from = c(estimate, std.error)) |> 
    rename(
        k_EF = estimate_EF, stderr_k_EF = std.error_EF, 
        k_PF = estimate_PF, stderr_k_PF = std.error_PF
    )

# Half-life in wide format
turn_t50_wd <- turn_t50 |> 
    pivot_wider(names_from = condition, values_from = c(Estimate, SE)) |> 
    rename(
        t50_EF = Estimate_EF, t50_PF = Estimate_PF, 
        stderr_t50_EF = SE_EF, stderr_t50_PF = SE_PF
    )

# Summary in wide format
turn_out <- turn_k_wd |> 
    left_join(turn_t50_wd) |> 
    mutate(log2FC = log2(k_EF) - log2(k_PF)) |> 
    left_join(
        turn_nested |> 
            mutate(n_peptide = map_int(data, \(x) length(unique(x[["peptide"]])))) |> 
            select(uniprot, site, n_peptide, pvalue)
    ) |> 
    mutate(adj_pvalue = p.adjust(pvalue, method = "BH")) |> 
    left_join(annot_turn) |> 
    select(
        uniprot, swissprot, gene_name, protein_name, site, n_peptide,
        k_PF, stderr_k_PF, k_EF, stderr_k_EF, 
        t50_PF, stderr_t50_PF, t50_EF, stderr_t50_EF, 
        log2FC, pvalue, adj_pvalue, localization
    ) |> 
    arrange(adj_pvalue, uniprot)

turn_sum <- turn_out |> 
    select(uniprot, gene_name, site, log2FC, pvalue, adj_pvalue, localization)

# Significance analysis of native protein turnover (Table S2B)
turn_out |> 
    filter(site == "unmodified") |> 
    select(-site) |> 
    arrange(adj_pvalue)

# Significance analysis of acetylated protein turnover (Table S5B)
turn_out |> 
    filter(site != "unmodified") |> 
    arrange(adj_pvalue)

# Plots for significance analysis -----------------------------------------

thresh_prac <- 1.5
thresh_stat <- 0.05

# LFQ of protein abundance
total_sum2 <- read_excel("data/ALD_MS.xlsx", sheet = "Table S1B", range = "A2:H1136") |> 
    mutate(
        loc = ifelse(
            localization %in% c("Mitochondrion", "MIM", "MM"),
            "Mitochondrion",
            ifelse(localization != "Cytoplasm", "Others", localization)
        ) |> factor(levels = c("Mitochondrion", "Cytoplasm", "Others"))
    )

# Volcano plot (Fig. 2B)
total_vlcn <- total_sum2 |> 
    ggplot(aes(log2FC, -log10(adj_pvalue))) + 
    geom_point(color = "darkgray", alpha = 0.7) +
    geom_point(aes(color = loc), alpha = 0.7, data = total_sum2 |> filter(adj_pvalue < thresh_stat)) +
    geom_hline(yintercept = -log10(0.05), color = "darkgray", linetype = "dashed") + 
    geom_vline(xintercept = c(1, -1) * log2(thresh_prac), color = "darkgray", linetype = "dashed") + 
    geom_text_repel(
        data = total_sum2 |> 
            filter(adj_pvalue < thresh_stat) |> 
            mutate(gene_name = case_match(
                gene_name, 
                "Hist1h2bp;H2bc12;H2bc4;Hist2h2bb;H2bc9;H2bc3;H2bc14;H2bc7" ~ "H2bc",
                "Rhoa;Rhoc" ~ "Rhoa/Rhoc", 
                "Timm8a1;Timm8a2" ~ "Timm8a1/2", 
                "Snx2;Snx1" ~ "Snx1/2",
                .default = gene_name
            )),
        aes(label = gene_name), max.overlaps = 25, show.legend = FALSE, size = 4
    ) + 
    scale_color_manual(values = cbPalette[c(6, 7, 4)]) +
    labs(
        x = bquote(~Log[2]~ "fold change in protein abundance (EF/PF)"),
        y = bquote(~-Log[10]~ "(adjusted p-value)"),
        color = ""
    ) + 
    coord_cartesian(xlim = c(-5, 4), ylim = c(0, 2.7))
total_vlcn + 
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.85), 
        legend.title = element_blank(),
        legend.background = element_rect(color = "darkgray")
    )

# LFQ of acetylation
aclevel_sum2 <- aclevel_sum |> 
    left_join(annot_ac |> select(uniprot, gene_name, localization)) |> 
    mutate(psite = str_c(gene_name, site, sep = "_")) |> 
    mutate(
        loc = ifelse(
            localization %in% c("Mitochondrion", "MIM", "Mitochondrial Matrix"), 
            "Mitochondrion", 
            ifelse(localization == "Cytoplasm", "Cytoplasm", "Others")
        ) |> factor(levels = c("Mitochondrion", "Cytoplasm", "Others"))
    )

# Volcano plot (Fig. S6A)
aclevel_sum2 |> 
    ggplot(aes(log2FC, -log10(adj_pvalue))) + 
    geom_point(color = "darkgray", alpha = 0.7) + 
    geom_point(aes(color = loc), alpha = 0.7, data = aclevel_sum2 |> filter(adj_pvalue < thresh_stat)) + 
    geom_hline(yintercept = -log10(0.05), color = "darkgray", linetype = "dashed") + 
    geom_vline(xintercept = c(1, -1) * log2(thresh_prac), color = "darkgray", linetype = "dashed") + 
    geom_text_repel(
        data = aclevel_sum2 |> filter(adj_pvalue < thresh_stat),
        aes(label = psite), max.overlaps = 25, show.legend = FALSE, size = 4
    ) + 
    scale_color_manual(values = cbPalette[c(6, 7, 4)]) + 
    labs(
        x = bquote(~Log[2]~ "fold change in acetylation (EF/PF)"),
        y = bquote(~-Log[10]~ "(adjusted p-value)"),
        color = ""
    )

# EtOH-induced change in turnover rate
turn_sum2 <- turn_sum |> 
    # mutate(modification = ifelse(site == "unmodified", "native", "acetylated")) |> 
    mutate(psite = ifelse(
        site == "unmodified", gene_name, str_c(gene_name, site, sep = "_")
    )) |> 
    mutate(
        loc = ifelse(
            str_detect(localization, "Mito") | localization == "MIM",
            "Mitochondrion", 
            ifelse(str_detect(localization, "Cyto"), "Cytoplasm", "Others")
        ) |> factor(levels = c("Mitochondrion", "Cytoplasm", "Others"))
    )

# Volcano plot for native protein turnover (Fig. 3D)
turn_sum_nt <- turn_sum2 |> filter(site == "unmodified")
turn_sum_nt |> 
    ggplot(aes(log2FC, -log10(adj_pvalue))) + 
    geom_point(color = "darkgray", alpha = 0.7) + 
    geom_point(
        aes(color = loc), alpha = 0.7, 
        data = turn_sum_nt |> filter(adj_pvalue < thresh_stat)
    ) + 
    geom_hline(yintercept = -log10(0.05), color = "darkgray", linetype = "dashed") + 
    geom_vline(xintercept = c(1, -1) * log2(thresh_prac), color = "darkgray", linetype = "dashed") + 
    geom_text_repel(
        data = turn_sum_nt |> 
            filter(adj_pvalue < thresh_stat, abs(log2FC) > log2(thresh_prac)),
        aes(label = psite), max.overlaps = 25, show.legend = FALSE, size = 4
    ) + 
    scale_color_manual(values = cbPalette[c(6, 7, 4)]) + 
    labs(
        x = bquote(~Log[2]~ "fold change in native protein turnover (EF/PF)"),
        y = bquote(~-Log[10]~ "(adjusted p-value)"),
        color = ""
    ) + 
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.85), 
        legend.title = element_blank(),
        legend.background = element_rect(color = "darkgray")
    )

# Volcano plot for acetylated protein turnover (Fig. 4D)
turn_sum_ac <- turn_sum2 |> filter(site != "unmodified")
turn_sum_ac |> 
    ggplot(aes(log2FC, -log10(adj_pvalue))) + 
    geom_point(color = "darkgray", alpha = 0.7) + 
    geom_point(aes(color = loc), alpha = 0.7, data = turn_sum_ac |> filter(adj_pvalue < thresh_stat)) + 
    geom_hline(yintercept = -log10(0.05), color = "darkgray", linetype = "dashed") + 
    geom_vline(xintercept = c(1, -1) * log2(thresh_prac), color = "darkgray", linetype = "dashed") + 
    geom_text_repel(
        data = turn_sum_ac |> filter(adj_pvalue < thresh_stat, abs(log2FC) > log2(thresh_prac)),
        aes(label = psite), max.overlaps = 40, show.legend = FALSE, size = 4
    ) + 
    scale_color_manual(values = cbPalette[c(6, 7, 4)]) + 
    labs(
        x = bquote(~Log[2]~ "fold change in acetylated protein turnover (EF/PF)"),
        y = bquote(~-Log[10]~ "(adjusted p-value)"),
        color = ""
    ) + 
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85), 
        legend.title = element_blank(),
        legend.background = element_rect(color = "darkgray")
    )

# Plots for protein-level turnover rate estimates -------------------------

k_box_cond <- turn_k |> 
    mutate(modification = ifelse(site == "unmodified", "native", "acetylated")) |> 
    group_by(uniprot, condition, modification) |> 
    summarise(k = mean(estimate), .groups = "drop") |> 
    mutate(modification = factor(modification, levels = c("native", "acetylated"))) |> 
    mutate(
        group = str_c(condition, modification) |> 
            fct_relevel("PFnative", "PFacetylated", "EFnative", "EFacetylated") |> 
            fct_recode("PF\nnative" = "PFnative", "PF\nacetylated" = "PFacetylated",
                       "EF\nnative" = "EFnative", "EF\nacetylated" = "EFacetylated")
    ) |> 
    left_join(annot_turn) |> 
    mutate(
        loc = ifelse(
            str_detect(localization, "Mito") | localization == "MIM",
            "Mitochondrion", 
            ifelse(str_detect(localization, "Cyto"), "Cytoplasm", "Others")
        ) |> factor(levels = c("Mitochondrion", "Cytoplasm", "Others"))
    )

# Protein-level turnover rate estimates (Fig. S1B)
set.seed(1)
k_box_cond |> 
    filter(modification == "native") |> 
    mutate(condition = factor(condition, levels = c("PF", "EF"))) |> 
    ggplot(aes(condition, k, color = condition)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(alpha = 0.6, width = 0.15, height = 0) + 
    scale_color_manual(values = pal_mod) + 
    labs(x = "", y = bquote("Turnover rate (day" ^-1*")")) + 
    theme_classic() + 
    theme(legend.position = "none")

# Protein-level turnover rate estimates in different compartments (Fig. S7B)
k_box_cond |> 
    filter(loc %in% c("Mitochondrion", "Cytoplasm")) |> 
    mutate(modification = factor(str_to_title(modification), levels = c("Native", "Acetylated"))) |>
    ggplot(aes(group, k)) + 
    geom_boxplot(aes(color = modification), outlier.shape = NA) + 
    geom_jitter(aes(color = modification), alpha = 0.6, width = 0.15, height = 0) + 
    scale_color_manual(values = pal_mod) + 
    labs(x = "", y = bquote("Turnover rate (day" ^-1*")")) + 
    facet_wrap(~ loc) + 
    theme(legend.position = "none")

# Turnover rates of native proteins in fatty acid oxidation (Fig. S2)
kk <- turn_k |> 
    filter(!str_detect(site, "K")) |> 
    select(-site) |> 
    left_join(annot_turn |> select(uniprot, gene_name)) |> 
    mutate(condition = factor(condition, c("PF", "EF")))
uu <- c("Q8VCW8", "Q8JZR0", "O35488", "Q91VA0", "Q8BMS1", "O08756", "Q61425", "Q8QZT1", 
        "Q8BWT1", "Q80XN0", "P52825", "P45952", "Q07417", "P42125", "P41216")
pvals <- turn_sum_nt |> 
    filter(uniprot %in% uu) |> 
    arrange(gene_name) |> 
    pull(pvalue) |> 
    annotate_pval() |> 
    str_replace("^p=0.[0-9]+", "")
kk |> filter(uniprot %in% uu) |> 
    ggplot(aes(x = gene_name, estimate, color = condition)) + 
    geom_col(position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_errorbar(
        aes(ymin = estimate-std.error, ymax = estimate+std.error), 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate("text", x = 1:length(uu), y = 0.022, label = pvals) + 
    scale_color_manual(values = pal_mod) +
    labs(x = "", y = bquote("Turnover rate ("~day^-1 ~")"), title = "Fatty acid oxidation", color = "") + 
    coord_cartesian(ylim = c(0, 0.023)) + 
    theme_classic()

# Turnover rates of native proteins in TCA cycle (Fig. S2)
uu <- c("Q99KI0", "Q9CZU6", "P97807", "P54071", "P08249", "Q60597")
pvals <- turn_sum_nt |> 
    filter(uniprot %in% uu) |> 
    arrange(gene_name) |> 
    pull(pvalue) |> 
    annotate_pval() |> 
    str_replace("^p=0.[0-9]+", "")
kk |> filter(uniprot %in% uu) |> 
    ggplot(aes(x = gene_name, estimate, color = condition)) + 
    geom_col(position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_errorbar(
        aes(ymin = estimate-std.error, ymax = estimate+std.error), 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate("text", x = 1:length(uu), y = 0.022, label = pvals) + 
    scale_color_manual(values = pal_mod) +
    labs(x = "", y = bquote("Turnover rate ("~day^-1 ~")"), title = "TCA cycle", color = "") + 
    coord_cartesian(ylim = c(0, 0.023)) + 
    theme_classic()

# Turnover rates of native proteins in ETC complex (Fig. S2)
uu <- c("Q03265", "P56480", "Q91VR2", "Q9DB20", "Q99LC3")
pvals <- turn_sum_nt |> 
    filter(uniprot %in% uu) |> 
    arrange(gene_name) |> 
    pull(pvalue) |> 
    annotate_pval() |> 
    str_replace("^p=0.[0-9]+", "")
kk |> filter(uniprot %in% uu) |> 
    ggplot(aes(x = gene_name, estimate, color = condition)) + 
    geom_col(position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_errorbar(
        aes(ymin = estimate-std.error, ymax = estimate+std.error), 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate("text", x = 1:length(uu), y = 0.022, label = pvals) + 
    scale_color_manual(values = pal_mod) +
    labs(x = "", y = bquote("Turnover rate ("~day^-1 ~")"), title = "ETC complex", color = "") + 
    coord_cartesian(ylim = c(0, 0.023)) + 
    theme_classic()

# Turnover rates of native proteins in alcohol metabolism (Fig. S2)
uu <- c("P00329", "P28474")
pvals <- turn_sum_nt |> 
    filter(uniprot %in% uu) |> 
    arrange(gene_name) |> 
    pull(pvalue) |> 
    annotate_pval() |> 
    str_replace("^p=0.[0-9]+", "")
kk |> filter(uniprot %in% uu) |> 
    ggplot(aes(x = gene_name, estimate, color = condition)) + 
    geom_col(position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_errorbar(
        aes(ymin = estimate-std.error, ymax = estimate+std.error), 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate("text", x = 1:length(uu), y = 0.022, label = pvals) + 
    scale_color_manual(values = pal_mod) +
    labs(x = "", y = bquote("Turnover rate ("~day^-1 ~")"), title = "Alcohol metabolism", color = "") + 
    coord_cartesian(ylim = c(0, 0.023)) + 
    theme_classic()

# Turnover rates of native proteins in urea cycle  (Fig. S2)
uu <- c("Q61176", "Q91YI0", "P16460", "Q8C196", "P11725", "P26443")
pvals <- turn_sum_nt |> 
    filter(uniprot %in% uu) |> 
    arrange(gene_name) |> 
    pull(pvalue) |> 
    annotate_pval() |> 
    str_replace("^p=0.[0-9]+", "")
kk |> filter(uniprot %in% uu) |> 
    ggplot(aes(x = gene_name, estimate, color = condition)) + 
    geom_col(position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_errorbar(
        aes(ymin = estimate-std.error, ymax = estimate+std.error), 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate("text", x = 1:length(uu), y = 0.022, label = pvals) + 
    scale_color_manual(values = pal_mod) +
    labs(x = "", y = bquote("Turnover rate ("~day^-1 ~")"), title = "Urea cycle", color = "") + 
    coord_cartesian(ylim = c(0, 0.023)) + 
    theme_classic()

# Turnover rates of native proteins in one-carbon metabolism (Fig. S2)
uu <- c("P24549", "P50431", "Q922D8", "P00375", "O35490", "Q91X83", "P50247",
        "Q8VCN5", "P97494")
pvals <- turn_sum_nt |> 
    filter(uniprot %in% uu) |> 
    arrange(gene_name) |> 
    pull(pvalue) |> 
    annotate_pval() |> 
    str_replace("^p=0.[0-9]+", "")
kk |> filter(uniprot %in% uu) |> 
    ggplot(aes(x = gene_name, estimate, color = condition)) + 
    geom_col(position = position_dodge(0.8), width = 0.7, fill = "white") +
    geom_errorbar(
        aes(ymin = estimate-std.error, ymax = estimate+std.error), 
        width = 0.2, position = position_dodge(0.8)
    ) + 
    annotate("text", x = 1:length(uu), y = 0.022, label = pvals) + 
    scale_color_manual(values = pal_mod) +
    labs(x = "", y = bquote("Turnover rate ("~day^-1 ~")"), title = "One-carbon metabolism", color = "") + 
    coord_cartesian(ylim = c(0, 0.023)) + 
    theme_classic()

# Plots for association between EtOH-induced changes ----------------------

# Association between native protein turnover and acetylated protein turnover (Fig. 4E)
cor_acturn_nturn <- turn_sum_ac |> 
    select(uniprot, site, psite, loc, log2FC_at = log2FC, adj_pvalue_at = adj_pvalue) |> 
    inner_join(
        turn_sum_nt |> select(uniprot, log2FC_nt = log2FC, adj_pvalue_nt = adj_pvalue)
    )
cor_acturn_nturn |> 
    ggplot(aes(log2FC_at, log2FC_nt)) + 
    annotate("rect", xmin = -1.2, xmax = 0, ymin = -1.2, ymax = 0, alpha = .05, fill = "red") + 
    geom_point(aes(color = loc), size = 2, alpha = 0.4) + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_text_repel(
        data = cor_acturn_nturn |> filter(adj_pvalue_at <= 0.05, adj_pvalue_nt <= 0.05),
        aes(color = loc, label = psite), size = 4, max.overlaps = 20, show.legend = FALSE
    ) +
    annotate(
        "text", x = -0.25, y = -1.05, size = 5, color = "darkred",
        label = "Decreased native & acetylated turnover",
        fontface = "italic"
    ) + 
    scale_color_manual(values = cbPalette[c(6, 7, 4)]) + 
    labs(
        x = bquote(~Log[2]~ "FC in acetylated protein turnover"),
        y = bquote(~Log[2]~ "FC in native protein turnover"),
        color = "",
    ) + 
    coord_cartesian(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) + 
    theme(legend.position = c(0.85, 0.25), legend.title = element_blank(),
          legend.background = element_rect(color = "darkgray"))

# cor_acturn_nturn |> distinct(uniprot)
# cor_acturn_nturn |> filter(log2FC_at < 0, log2FC_nt < 0) |> distinct(uniprot, loc) |> count(loc)
# cor_acturn_nturn |> filter(adj_pvalue_at < 0.05, adj_pvalue_nt < 0.05) 

# Association between acetylation level and acetylated protein turnover (Fig. S9A)
cor_acturn_aclevel <- turn_sum_ac |>
    inner_join(
        aclevel_sum |> select(uniprot, site, log2FC_gbl = log2FC, adj_pvalue_gbl = adj_pvalue)
    )
cor_acturn_aclevel |> 
    ggplot(aes(log2FC, log2FC_gbl, color = loc)) + 
    geom_point(size = 2, alpha = 0.4) + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_text_repel(
        data = cor_acturn_aclevel |> filter(adj_pvalue <= 0.05, adj_pvalue_gbl <= 0.05),
        aes(label = psite), size = 3.5, max.overlaps = 20, show.legend = FALSE
    ) +
    scale_color_manual(values = cbPalette[c(6, 7, 4)]) + 
    labs(
        x = bquote(~Log[2]~ "FC in acetylated protein turnover"),
        y = bquote(~Log[2]~ "FC in acetylation"),
        color = ""
    )

# Association between acetylated protein turnover and protein abundance (Fig. S9B)
cor_acturn_total <- turn_sum_ac |>
    select(uniprot, site, psite, loc, log2FC_at = log2FC, adj_pvalue_at = adj_pvalue) |> 
    inner_join(
        total_sum |> select(uniprot, log2FC_total = log2FC, adj_pvalue_total = adj_pvalue)
    )
cor_acturn_total |> 
    ggplot(aes(log2FC_at, log2FC_total, color = loc)) + 
    geom_point(size = 2, alpha = 0.4) + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_text_repel(
        data = cor_acturn_total |> filter(adj_pvalue_at <= 0.05, adj_pvalue_total <= 0.05),
        aes(label = psite), size = 3.5, max.overlaps = 20, show.legend = FALSE
    ) +
    scale_color_manual(values = cbPalette[c(6, 7, 4)]) + 
    labs(
        x = bquote(~Log[2]~ "FC in acetylated protein turnover"),
        y = bquote(~Log[2]~ "FC in protein abundance"),
        color = ""
    )

# Association between acetylation level and protein abundance (Fig. S6B)
cor_aclevel_total <- aclevel_sum2 |> 
    select(uniprot, site, psite, loc, log2FC_aclevel = log2FC, adj_pvalue_aclevel = adj_pvalue) |> 
    inner_join(
        total_sum |> select(uniprot, log2FC_total = log2FC, adj_pvalue_total = adj_pvalue)
    )
cor_aclevel_total |> 
    ggplot(aes(log2FC_total, log2FC_aclevel, color = loc)) + 
    geom_point(size = 2, alpha = 0.4) + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_text_repel(
        data = cor_aclevel_total |> filter(adj_pvalue_aclevel <= 0.05, adj_pvalue_total <= 0.05),
        aes(label = psite), size = 3.5, max.overlaps = 20, show.legend = FALSE
    ) +
    scale_color_manual(values = cbPalette[c(6, 7, 4)]) + 
    labs(
        x = bquote(~Log[2]~ "FC in protein abundance"),
        y = bquote(~Log[2]~ "FC in acetylation"),
        color = ""
    )

# Effect of acetylation on total protein quantification? (Fig. S6C)
total_long |> 
    group_by(uniprot, gene_name, condition) |> 
    summarise(log2inty = mean(log2inty, na.rm = TRUE), .groups = "drop") |> 
    mutate(acetylated = ifelse(uniprot %in% aclevel_sum$uniprot, "yes", "no")) |> 
    mutate(
        condition = ifelse(condition == "EF", "Ethanol-fed", "Pair-fed") |> 
            factor(levels = c("Pair-fed", "Ethanol-fed"))
    ) |> 
    ggplot(aes(acetylated, log2inty, color = condition)) + 
    geom_boxplot() + 
    scale_color_manual(values = c("darkgray", "darkred")) +
    labs(x = "Acetylation detected and quantified", 
         y = bquote(~Log[2]~ "LFQ of protein abundance")) + 
    facet_wrap(~ condition)

# Heatmap (Fig. S9C)
library(gplots)
hmcol <- bluered(100)

log2FC_all <- cor_acturn_nturn |> 
    inner_join(cor_aclevel_total |> select(psite, contains("log2FC"), contains("adj"))) |> 
    select(psite, loc, contains("log2FC")) |> 
    arrange(loc)

log2FC_all |> count(loc)

log2FC_mat <- as.matrix(log2FC_all[, 3:6])
rownames(log2FC_mat) <- log2FC_all$psite
colnames(log2FC_mat) <- c("Acetylated protein\nturnover", "Native protein\nturnover", "Acetylation", "Protein\nabundance")

heatmap.2(
    x = t(log2FC_mat), 
    dendrogram = "col", 
    col = hmcol, ColSideColors = rep(cbPalette[c(6, 7, 4)], c(111, 13, 4)), 
    lmat=rbind(c(5,4,0), c(0,1,0), c(3,2,0)), lwid = c(1.0, 6, 0.6), lhei = c(1.5, 0.2, 5),
    trace = "none", scale = "none", density.info = "none", cexRow = 1.5, cexCol = 0.6
)

# Functional enrichment analysis ------------------------------------------

pat_ac <- paste0(
    "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]",
    "([A-Z][A-Z0-9]{2}[0-9]){1,2})"
)

# Background 1 (from PMID: 34836951, MassIVE ID: MSV000086426)
atlas_liver <- read_tsv(
    "data/misc/MASSIVE-COMPLETE-5abc0707-display_quant_results-main.tsv",
    col_select = contains("Protein")
) |> 
    rename(protein_group = `Protein Groups`) |> 
    mutate(uniprots = str_extract_all(protein_group, pat_ac))
bg_liver <- unique(unlist(atlas_liver$uniprots))

# Background 2 (from PMID: 36058520)
bg_liver_total <- read_csv(
    "data/misc/PMID36058520-mmc3.csv",
    col_select = c(1, contains("LIV"))
) |> 
    filter(!is.na(LIV1) | !is.na(LIV2) | !is.na(LIV3) | !is.na(LIV4)) |> 
    pull(1)

bg_liver_full <- unique(c(bg_liver_total, bg_liver))

# Proteins with significant changes
sig_nt <- turn_sum_nt |> 
    filter(adj_pvalue < 0.05) |> 
    select(uniprot, log2FC)

sig_at <- turn_sum_ac |> 
    filter(adj_pvalue < 0.05) |> 
    group_by(uniprot) |> 
    summarise(log2FC = mean(log2FC))

# FC vectors
fc_nt <- sig_nt$log2FC
fc_at <- sig_at$log2FC

names(fc_nt) <- sig_nt$uniprot
names(fc_at) <- sig_at$uniprot

# compare differential protein lists
dp_list <- list(
    increased_native_turnover = names(fc_nt)[fc_nt > 0],
    decreased_native_turnover = names(fc_nt)[fc_nt < 0],
    increased_acetylated_turnover = names(fc_at)[fc_at > 0],
    decreased_acetylated_turnover = names(fc_at)[fc_at < 0]
)
dp_vec <- unlist(dp_list)

dp_df <- data.frame(protein = unname(dp_vec), group = names(dp_vec)) |> 
    mutate(
        char = ifelse(str_detect(group, "native"), "Native turnover", "Acetylated turnover"),
        change = ifelse(str_detect(group, "increased"), "Increased", "Decreased")
    )

# Load necessary packages
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

# Functional enrichment analysis based on proteins with significant changes in 
# turnover rate
cgo_nt <- compareCluster(
    protein ~ char + change, data = dp_df |> filter(char == "Native turnover"), 
    fun = "enrichGO", 
    universe = bg_liver_full, 
    keyType = "UNIPROT", 
    OrgDb = org.Mm.eg.db, 
    ont = "BP", 
    pAdjustMethod = "BH", 
    qvalueCutoff = 0.05,
    readable = TRUE
)

# Functional enrichment analysis based on proteins whose acetylated forms have 
# significant changes in turnover rate
cgo_at <- compareCluster(
    protein ~ char + change, data = dp_df |> filter(char == "Acetylated turnover"), 
    fun = "enrichGO", 
    universe = bg_liver_full, 
    keyType = "UNIPROT", 
    OrgDb = org.Mm.eg.db, 
    ont = "BP", 
    pAdjustMethod = "BH", 
    qvalueCutoff = 0.05,
    readable = TRUE
)

# Remove redundant terms
cgo_nt2 <- simplify(cgo_nt)
cgo_at2 <- simplify(cgo_at)

# Significant GO terms (Table S2C)
as.data.frame(cgo_nt2) |> dplyr::select(-c(char, -change))

# Significant GO terms (Table S5C)
as.data.frame(cgo_at2) |> dplyr::select(-c(char, -change))

# Dot plot (Fig. 3E)
dotplot(cgo_nt2, x = "change", showCategory = 6) + 
    facet_grid(~ char) + 
    scale_color_gradientn(
        colours=c("#b3eebe", "#46bac2", "#371ea3"), guide=guide_colorbar(reverse=TRUE, order=1)
    ) +
    labs(x = "", size = "Ratio")

# Dot plot (Fig. 4F)
dotplot(cgo_at2, x = "change") + 
    facet_grid(~ char) + 
    scale_color_gradientn(
        colours=c("#b3eebe", "#46bac2", "#371ea3"), guide=guide_colorbar(reverse=TRUE, order=1)
    ) +
    labs(x = "", size = "Ratio")

# Pooled analyses
go_nt <- enrichGO(
    names(fc_nt),
    universe = bg_liver_full,
    keyType = "UNIPROT",
    OrgDb = org.Mm.eg.db,
    ont="BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
)

go_at <- enrichGO(
    names(fc_at),
    universe = bg_liver_full,
    keyType = "UNIPROT",
    OrgDb = org.Mm.eg.db,
    ont="BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
)

# Remove redundant terms
go_nt2 <- simplify(go_nt)
go_at2 <- simplify(go_at)

# Pairwise
go_nt2 <- pairwise_termsim(go_nt2)
go_at2 <- pairwise_termsim(go_at2)

# Gene-concept network (Fig. 3F)
cnetplot(go_nt2, color.params = list(foldChange = fc_nt), showCategory = 8)

# Gene-concept network (Fig. 4G)
cnetplot(go_at2, color.params = list(foldChange = fc_at), showCategory = 8)
