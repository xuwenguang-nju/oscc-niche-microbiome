#!/usr/bin/env Rscript
# A3: overall survival analysis for the intratumoral microbiome in TCGA-HNSC.
#
# Uses the `survival` package (3.5.8, on the host) so log-rank and Cox results are standard.
#
# WHY STAGE ENTERS AS AN ORDINAL, NOT A FACTOR
# Stage I has 7 patients and 0 deaths, so a factor coding hits complete separation: every
# other level's HR against Stage I goes to ~1e7 with an infinite CI and p=0.997. That is a
# numerical artefact, not a finding. Coding stage as a continuous ordinal (1-4) removes the
# zero-event reference level, uses all 106 patients, and is interpreted as HR per one-stage
# increase. A binary early(I-II)/late(III-IV) split is fitted as a sensitivity check.
#
# TWO QUESTIONS, KEPT SEPARATE
#   1. prognostic?  -> continuous Cox per SD, plus KM by median for the figure.
#   2. independent of stage/age/sex? -> multivariable Cox + likelihood-ratio test of whether
#      adding the microbiome term improves on the clinical-only model. Without the LRT a
#      "significant" adjusted HR is easy to over-read.
#
# Samples are collapsed to one row per patient before modelling.
#
# Run: Rscript s11_survival.R
suppressPackageStartupMessages({
  library(survival)
  library(jsonlite)
})

d <- read.csv("master.csv", stringsAsFactors = FALSE)
cat("samples:", nrow(d), "\n")

metrics <- c("shannon", "richness", "parvimonas_rel")
first_non_na <- function(x) {
  v <- x[!is.na(x)]
  if (length(v)) v[1] else NA
}

agg <- aggregate(d[, metrics], by = list(case_id = d$case_id), FUN = mean, na.rm = TRUE)
clin_cols <- c("os_time", "os_event", "age", "sex", "stage", "grade_num",
               "site_group", "smoker", "hpv")
first <- aggregate(d[, clin_cols], by = list(case_id = d$case_id), FUN = first_non_na)
p <- merge(agg, first, by = "case_id")
p <- p[!is.na(p$os_time) & p$os_time > 0, ]

# stage: ordinal + binary
p$stage_num <- as.numeric(p$stage)
p$stage_bin <- ifelse(p$stage_num <= 2, "early (I-II)",
                      ifelse(p$stage_num >= 3, "late (III-IV)", NA))
cat("patients with OS:", nrow(p), "| deaths:", sum(p$os_event),
    "| with stage:", sum(!is.na(p$stage_num)), "\n")
cat("\ndeaths by stage:\n")
print(table(stage = p$stage_num, death = p$os_event, useNA = "ifany"))

for (m in metrics) p[[paste0(m, "_z")]] <- as.numeric(scale(p[[m]]))

# ---- 1. univariable Cox, per SD -------------------------------------------
uni <- list()
for (m in metrics) {
  z <- paste0(m, "_z")
  fit <- coxph(as.formula(paste0("Surv(os_time, os_event) ~ ", z)), data = p)
  s <- summary(fit)
  uni[[m]] <- data.frame(
    metric = m, HR_per_SD = round(s$coefficients[1, "exp(coef)"], 3),
    lower95 = round(s$conf.int[1, "lower .95"], 3),
    upper95 = round(s$conf.int[1, "upper .95"], 3),
    p = signif(s$coefficients[1, "Pr(>|z|)"], 4),
    n = fit$n, events = fit$nevent, concordance = round(s$concordance[1], 3))
}
U <- do.call(rbind, uni)

# ---- 2. KM by median split ------------------------------------------------
km <- list()
for (m in metrics) {
  x <- p[[m]]
  if (sum(!is.na(x)) < 20) next
  p$.grp <- ifelse(x > median(x, na.rm = TRUE), "high", "low")
  sd <- survdiff(Surv(os_time, os_event) ~ .grp, data = p)
  pv <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  km[[m]] <- data.frame(metric = m, n_high = sum(p$.grp == "high"),
                        n_low = sum(p$.grp == "low"),
                        chisq = round(sd$chisq, 3), p_logrank = signif(pv, 4))
}
p$.grp <- NULL
K <- do.call(rbind, km)

cat("\n=== 1. univariable Cox (continuous, per 1 SD) ===\n"); print(U, row.names = FALSE)
cat("\n=== 2. Kaplan-Meier by median split (log-rank) ===\n"); print(K, row.names = FALSE)

# ---- 3. clinical reference model (stage ordinal) --------------------------
dc <- p[!is.na(p$stage_num) & !is.na(p$age) & !is.na(p$sex), ]
fit_clin <- coxph(Surv(os_time, os_event) ~ stage_num + age + sex, data = dc)
sc <- summary(fit_clin)
cat("\n=== 3. clinical reference model (stage ordinal + age + sex) ===\n")
print(data.frame(term = rownames(sc$coefficients),
                 HR = round(sc$coefficients[, "exp(coef)"], 3),
                 lower95 = round(sc$conf.int[, "lower .95"], 3),
                 upper95 = round(sc$conf.int[, "upper .95"], 3),
                 p = signif(sc$coefficients[, "Pr(>|z|)"], 4)), row.names = FALSE)
cat("concordance:", round(sc$concordance[1], 3), " n:", fit_clin$n,
    " events:", fit_clin$nevent, "\n")

# ---- 4. does the microbiome add anything? (LRT vs clinical-only) ----------
mv <- list(); lrt <- list()
for (m in metrics) {
  z <- paste0(m, "_z")
  dd <- p[!is.na(p$stage_num) & !is.na(p$age) & !is.na(p$sex), ]
  f1 <- as.formula(paste0("Surv(os_time, os_event) ~ ", z, " + stage_num + age + sex"))
  fit1 <- coxph(f1, data = dd)
  fit0 <- coxph(Surv(os_time, os_event) ~ stage_num + age + sex, data = dd)
  s1 <- summary(fit1)
  mv[[m]] <- data.frame(
    metric = m, term = rownames(s1$coefficients),
    HR = round(s1$coefficients[, "exp(coef)"], 3),
    lower95 = round(s1$conf.int[, "lower .95"], 3),
    upper95 = round(s1$conf.int[, "upper .95"], 3),
    p = signif(s1$coefficients[, "Pr(>|z|)"], 4), row.names = NULL)
  lrt[[m]] <- data.frame(
    metric = m, loglik_clin = round(fit0$loglik[2], 3),
    loglik_clin_plus_metric = round(fit1$loglik[2], 3),
    LRT_chisq = round(2 * (fit1$loglik[2] - fit0$loglik[2]), 3),
    LRT_p = signif(pchisq(2 * (fit1$loglik[2] - fit0$loglik[2]), df = 1,
                          lower.tail = FALSE), 4),
    AIC_clin = round(AIC(fit0), 2), AIC_clin_plus_metric = round(AIC(fit1), 2),
    n = fit1$n, events = fit1$nevent)
}
MV <- do.call(rbind, mv); LRT <- do.call(rbind, lrt)
cat("\n=== 4. multivariable Cox (adjusted for stage ordinal, age, sex) ===\n")
print(MV, row.names = FALSE)
cat("\n=== 4b. likelihood-ratio test: clinical-only vs clinical + microbiome ===\n")
print(LRT, row.names = FALSE)

# ---- 5. sensitivity: binary stage ----------------------------------------
ds <- p[!is.na(p$stage_bin) & !is.na(p$age) & !is.na(p$sex), ]
ds$stage_bin <- factor(ds$stage_bin)
sens <- list()
for (m in c("shannon", "richness")) {
  z <- paste0(m, "_z")
  fit <- coxph(as.formula(paste0("Surv(os_time, os_event) ~ ", z,
                                 " + stage_bin + age + sex")), data = ds)
  s <- summary(fit)
  sens[[m]] <- data.frame(metric = m, term = rownames(s$coefficients),
                          HR = round(s$coefficients[, "exp(coef)"], 3),
                          p = signif(s$coefficients[, "Pr(>|z|)"], 4), row.names = NULL)
}
S <- do.call(rbind, sens)
cat("\n=== 5. sensitivity: stage as early(I-II) vs late(III-IV) ===\n")
print(S, row.names = FALSE)

write.csv(U, "s11_survival_univariable.csv", row.names = FALSE)
write.csv(K, "s11_survival_km.csv", row.names = FALSE)
write.csv(MV, "s11_survival_multivariable.csv", row.names = FALSE)
write.csv(LRT, "s11_survival_lrt.csv", row.names = FALSE)
write.csv(S, "s11_survival_stage_sensitivity.csv", row.names = FALSE)

kmfit <- survfit(Surv(os_time, os_event) ~ grp,
                 data = transform(p, grp = ifelse(
                   shannon > median(shannon, na.rm = TRUE),
                   "high shannon", "low shannon")))
saveRDS(kmfit, "s11_km_shannon.rds")

write(toJSON(list(n_samples = nrow(d), n_patients = nrow(p),
                  events = as.integer(sum(p$os_event)),
                  n_stage = sum(!is.na(p$stage_num)),
                  stage_coding = "ordinal 1-4 (factor coding hit complete separation at stage I)",
                  metrics = metrics), auto_unbox = TRUE, pretty = TRUE),
      "s11_survival_status.json")
cat("\nsaved: s11_survival_univariable / _km / _multivariable / _lrt / _stage_sensitivity.csv\n")
