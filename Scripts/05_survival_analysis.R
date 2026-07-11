## 05_survival_analysis.R
## Build the survival dataset, run Kaplan-Meier and Cox regression for the
## complement-evasion score (unadjusted and stage-adjusted), validate the immune
## score against survival, and test its correlation with complement-evasion score.

library(survival)
library(survminer)
library(dplyr)

## --- Build survival time/event variables from TCGA clinical data ---
clinical <- colData(tumor_only)

clinical$time <- ifelse(clinical$vital_status == "Dead",
                         clinical$days_to_death,
                         clinical$days_to_last_follow_up)
clinical$event <- ifelse(clinical$vital_status == "Dead", 1, 0)

surv_df <- data.frame(
  time = clinical$time,
  event = clinical$event,
  complement_score = tcga_complement_score
)
surv_df <- surv_df[!is.na(surv_df$time) & surv_df$time > 0, ]   # n = 518

## --- Kaplan-Meier: complement-evasion score (median split) ---
surv_df$group <- ifelse(surv_df$complement_score > median(surv_df$complement_score),
                         "High", "Low")
fit <- survfit(Surv(time, event) ~ group, data = surv_df)

km_plot <- ggsurvplot(fit, data = surv_df, pval = TRUE, risk.table = TRUE,
                       title = "Overall Survival by Complement-Evasion Score",
                       legend.title = "Group",
                       palette = c("#E7B800", "#2E9FDF"))

pdf("KM_plot_aligned.pdf", width = 8, height = 7)
print(km_plot)   # print.ggsurvplot handles plot/risk-table alignment correctly
dev.off()
# Figure 2 in manuscript

## --- Cox regression: complement-evasion score, unadjusted ---
cox_fit <- coxph(Surv(time, event) ~ complement_score, data = surv_df)
summary(cox_fit)

## --- Stage-adjusted Cox model. NOTE: modeling all 7 AJCC substages separately
## produced non-convergent, numerically unstable coefficients (sparse substage
## counts); stage was collapsed to Early (I-II) vs. Late (III-IVC) to fix this. ---
surv_df$stage <- clinical$ajcc_pathologic_stage[match(rownames(surv_df), rownames(clinical))]
surv_df$stage_group <- case_when(
  surv_df$stage %in% c("Stage I", "Stage II") ~ "Early",
  surv_df$stage %in% c("Stage III", "Stage IVA", "Stage IVB", "Stage IVC") ~ "Late",
  TRUE ~ NA_character_
)
table(surv_df$stage_group, useNA = "always")

cox_fit_adjusted2 <- coxph(Surv(time, event) ~ complement_score + stage_group, data = surv_df)
summary(cox_fit_adjusted2)

## --- Immune score: positive control (should predict survival) + correlation with
## complement-evasion score (Section 3.4, Table 1) ---
surv_df$immune_score <- tcga_immune_score[match(rownames(surv_df), colnames(tumor_only))]

cor.test(surv_df$complement_score, surv_df$immune_score, method = "spearman")

cox_immune_check <- coxph(Surv(time, event) ~ immune_score + stage_group, data = surv_df)
summary(cox_immune_check)   # validates immune score: significant survival association expected
