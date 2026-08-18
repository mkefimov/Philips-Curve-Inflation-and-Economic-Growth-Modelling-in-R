pkgs <- c("readxl", "dplyr", "ggplot2", "nlme", "sandwich", "lmtest",
          "car", "gmm", "tseries", "forecast", "patchwork", "ggrepel", "tidyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(readxl); library(dplyr);    library(ggplot2); library(nlme)
library(sandwich); library(lmtest); library(car);     library(gmm)
library(tseries);  library(forecast); library(patchwork); library(ggrepel)
library(tidyr)


theme_kp <- theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


raw <- read_excel("Phillips_Curve_Russia_1991_2026.xlsx",
                  sheet = "Phillips Curve Data")

df <- raw %>%
  rename(
    year        = 1,  cpi        = 2,  output_gap = 6,
    inf_expect  = 7,  key_rate   = 8,  reer       = 11,
    budget_def  = 12,                  
    cpi_lag1    = 13, d_reer     = 14
  ) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  filter(!is.na(year), year >= 1991) %>%
  arrange(year)

df <- df %>%
  mutate(
    d_sanctions   = as.integer(year %in% c(2014, 2022, 2023)),
    d_covid       = as.integer(year == 2020),
    d_budget      = budget_def - lag(budget_def, 1),    
    fiscal_imp    = -d_budget,
        fiscal_tight  = as.integer(!is.na(d_budget) & d_budget > 0)
  )


df_model <- df %>%
  filter(!is.na(cpi) & !is.na(cpi_lag1) & !is.na(output_gap) &
           !is.na(inf_expect) & !is.na(reer) & !is.na(d_reer) &
           !is.na(key_rate)  & !is.na(budget_def)) %>%
  mutate(year = as.integer(year))

cat("Размер выборки для моделирования:", nrow(df_model),
    "наблюдений, годы:", min(df_model$year), "–", max(df_model$year), "\n\n")

cat("\n", strrep("=", 70), "\n")
cat("Тесты на стационарность (ADF и KPSS)\n")
cat(strrep("=", 70), "\n\n")

unit_root_vars <- c("cpi", "cpi_lag1", "inf_expect", "output_gap",
                    "d_reer", "key_rate", "reer", "budget_def", "fiscal_imp")

run_ur_tests <- function(varname, data) {
  x <- na.omit(data[[varname]])
  adf  <- tryCatch(
    suppressWarnings(adf.test(x, k = trunc((length(x) - 1)^(1/3)))),
    error = function(e) list(statistic = NA, p.value = NA)
  )
  kpss <- tryCatch(
    suppressWarnings(kpss.test(x, null = "Level")),
    error = function(e) list(statistic = NA, p.value = NA)
  )
  data.frame(
    Переменная   = varname,
    N            = length(x),
    ADF_stat     = round(adf$statistic,  3),
    ADF_p        = round(adf$p.value,    4),
    ADF_вывод    = ifelse(!is.na(adf$p.value) & adf$p.value < 0.05,
                          "стационарен", "нестационарен"),
    KPSS_stat    = round(kpss$statistic, 3),
    KPSS_p       = round(kpss$p.value,   4),
    KPSS_вывод   = ifelse(!is.na(kpss$p.value) & kpss$p.value > 0.05,
                          "стационарен", "нестационарен"),
    check.names  = FALSE
  )
}

ur_table <- lapply(unit_root_vars, run_ur_tests, data = df_model) %>% bind_rows()

cat("Результаты тестов на единичный корень:\n")
cat("ADF  — H₀: ряд содержит единичный корень (нестационарен)\n")
cat("KPSS — H₀: ряд стационарен\n\n")
print(ur_table[ , c("Переменная","N","ADF_stat","ADF_p","ADF_вывод",
                    "KPSS_stat","KPSS_p","KPSS_вывод")], row.names = FALSE)

nonstat_vars <- ur_table %>%
  filter(ADF_вывод == "нестационарен" | KPSS_вывод == "нестационарен") %>%
  pull(Переменная)

if (length(nonstat_vars) == 0) {
  cat("\n Все переменные стационарны по обоим тестам.\n")
} else {
  cat(sprintf("\n Потенциально нестационарные переменные: %s\n",
              paste(nonstat_vars, collapse = ", ")))
  cat("  В данной работе применяются HAC-ошибки и GLS AR(1).\n")
  cat("  Для budget_def при нестационарности используется fiscal_imp \n")
}

ur_plot_df <- df_model %>%
  select(year, cpi, inf_expect, output_gap, budget_def, fiscal_imp) %>%
  pivot_longer(-year, names_to = "var", values_to = "val") %>%
  mutate(var = dplyr::recode(var,
                             cpi        = "ИПЦ (cpi)",
                             inf_expect = "Инф. ожидания",
                             output_gap = "Разрыв выпуска",
                             budget_def = "Бюджетный баланс (% ВВП)",
                             fiscal_imp = "Фискальный импульс (delta баланса)"
  ))

p_ur <- ggplot(ur_plot_df, aes(x = year, y = val)) +
  geom_line(color = "#4575b4", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  facet_wrap(~var, scales = "free_y", ncol = 2) +
  labs(
    title    = "Динамика ключевых переменных модели",
    subtitle = "Включая бюджетный баланс и фискальный импульс Минфина РФ",
    x = "Год", y = NULL
  ) + theme_kp
print(p_ur)
ggsave("fig0_stationarity_dynamics.png", p_ur,
       width = 11, height = 8, dpi = 150, bg = "white")

cat("\n", strrep("=", 70), "\n")
cat("Базовая гибридная кривая Филлипса (HAC-ошибки)\n")
cat(strrep("=", 70), "\n\n")

ols_base <- lm(cpi ~ inf_expect + cpi_lag1 + output_gap, data = df_model)
hac_base <- coeftest(ols_base, vcov = NeweyWest(ols_base, lag = 2, prewhite = FALSE))

cat("OLS (обычные ошибки)\n"); print(summary(ols_base))
cat("\nHAC-робастные ошибки (Newey-West, lag=2)\n"); print(hac_base)

coef_b <- coef(ols_base)
cat(sprintf("\nR^2 = %.4f | γ_f = %.4f | γ_b = %.4f | λ = %.4f | γ_f+γ_b = %.4f\n",
            summary(ols_base)$r.squared,
            coef_b["inf_expect"], coef_b["cpi_lag1"], coef_b["output_gap"],
            coef_b["inf_expect"] + coef_b["cpi_lag1"]))

cat("\n", strrep("=", 70), "\n")
cat("Поэтапное добавление переменных (включая бюджет Минфина)\n")
cat(strrep("=", 70), "\n\n")

models <- list(
  "M1: Базовая"          = lm(cpi ~ inf_expect + cpi_lag1 + output_gap,
                              data = df_model),
  "M2: + d_reer"         = lm(cpi ~ inf_expect + cpi_lag1 + output_gap + d_reer,
                              data = df_model),
  "M3: + key_rate"       = lm(cpi ~ inf_expect + cpi_lag1 + output_gap + d_reer + key_rate,
                              data = df_model),
  "M4: + d_sanctions"    = lm(cpi ~ inf_expect + cpi_lag1 + output_gap + d_reer +
                                key_rate + d_sanctions,
                              data = df_model),
  "M5: + d_covid"        = lm(cpi ~ inf_expect + cpi_lag1 + output_gap + d_reer +
                                key_rate + d_sanctions + d_covid,
                              data = df_model),
  "M6: + budget_def"     = lm(cpi ~ inf_expect + cpi_lag1 + output_gap + d_reer +
                                key_rate + d_sanctions + d_covid + budget_def,
                              data = df_model),
  "M7: + fiscal_imp"     = lm(cpi ~ inf_expect + cpi_lag1 + output_gap + d_reer +
                                key_rate + d_sanctions + d_covid + fiscal_imp,
                              data = filter(df_model, !is.na(fiscal_imp))),
  "M8: полная фискальная"= lm(cpi ~ inf_expect + cpi_lag1 + output_gap + d_reer +
                                key_rate + d_sanctions + d_covid +
                                budget_def + fiscal_imp,
                              data = filter(df_model, !is.na(fiscal_imp)))
)

step_tbl <- lapply(names(models), function(nm) {
  m  <- models[[nm]]
  hc <- coeftest(m, vcov = NeweyWest(m, lag = 2, prewhite = FALSE))
  vf <- if (length(coef(m)) > 2) tryCatch(vif(m), error = function(e) rep(NA, length(coef(m))-1)) else NA
  max_vif <- if (!all(is.na(vf))) max(vf, na.rm = TRUE) else NA
  data.frame(
    Модель   = nm, N = nobs(m),
    R2       = round(summary(m)$r.squared,     4),
    adj_R2   = round(summary(m)$adj.r.squared, 4),
    AIC      = round(AIC(m), 2),
    BIC      = round(BIC(m), 2),
    max_VIF  = round(max_vif, 2)
  )
}) %>% bind_rows()

cat("Сводная таблица пошагового добавления переменных (M1–M8):\n")
print(step_tbl, row.names = FALSE)

cat("\nHAC-коэффициенты для фискальных спецификаций (M6, M7, M8):\n")
for (nm in c("M6: + budget_def", "M7: + fiscal_imp", "M8: полная фискальная")) {
  cat(sprintf("\n  >> %s\n", nm))
  print(coeftest(models[[nm]], vcov = NeweyWest(models[[nm]], lag = 2, prewhite = FALSE)))
}

cat("\nVIF для полной фискальной модели M8:\n")
print(round(vif(models[["M8: полная фискальная"]]), 3))

cat("\n", strrep("=", 70), "\n")
cat("Оценка фискальной политики Минфина РФ\n")
cat(strrep("=", 70), "\n\n")

m_fiscal <- models[["M6: + budget_def"]]
hac_fiscal <- coeftest(m_fiscal, vcov = NeweyWest(m_fiscal, lag = 2, prewhite = FALSE))

cat("Модель M6 (базовая + бюджетный баланс)\n")
print(hac_fiscal)

coef_f   <- hac_fiscal["budget_def", "Estimate"]
se_f     <- hac_fiscal["budget_def", "Std. Error"]
pval_f   <- hac_fiscal["budget_def", "Pr(>|t|)"]
sign_str <- ifelse(pval_f < 0.05, "Значим (p < 0.05)", "Незначим (p ≥ 0.05)")

cat(sprintf("\n──────────────────────────────────────────────────────────────\n"))
cat(sprintf("  Коэффициент budget_def: %.4f (SE = %.4f, p = %.4f) — %s\n",
            coef_f, se_f, pval_f, sign_str))
cat(sprintf("──────────────────────────────────────────────────────────────\n\n"))

cat("Интерпретация\n")
if (coef_f < 0) {
  cat(sprintf(
    "Знак отрицательный (%.4f): рост бюджетного профицита на 1 п.п. ВВП\n",
    coef_f))
  cat(sprintf("ассоциируется со снижением инфляции на %.4f п.п.\n", abs(coef_f)))
  cat("Интерпретация: ужесточение бюджетной политики (сокращение дефицита)\n")
  cat("оказывает дезинфляционный эффект — подтверждает антиинфляционную роль\n")
  cat("Минфина через бюджетную консолидацию.\n")
} else {
  cat(sprintf(
    "Знак положительный (%.4f): рост дефицита на 1 п.п. ВВП\n", coef_f))
  cat(sprintf("ассоциируется с ростом инфляции на %.4f п.п.\n", abs(coef_f)))
  cat("Интерпретация: проциклическая или стимулирующая бюджетная политика\n")
  cat("оказывает инфляционное давление. Минфину необходима консолидация.\n")
}

cat("Хронология бюджетной политики Минфина РФ (по данным модели):\n\n")
fiscal_history <- df_model %>%
  select(year, budget_def, d_budget, fiscal_tight) %>%
  mutate(
    Политика = case_when(
      !is.na(fiscal_tight) & fiscal_tight == 1 ~ "Ужесточение (консолидация)",
      !is.na(d_budget) & d_budget < -1          ~ "Значимое смягчение",
      !is.na(d_budget) & d_budget <  0          ~ "Умеренное смягчение",
      TRUE                                       ~ "Нейтральная/профицит"
    )
  ) %>%
  select(Год = year, `Баланс (% ВВП)` = budget_def,
         `delta баланса` = d_budget, Политика)
print(fiscal_history, row.names = FALSE)

cat("\nСредняя инфляция по режимам бюджетной политики:\n")
df_model %>%
  mutate(
    Режим = case_when(
      budget_def > 0 ~ "Профицит",
      budget_def >= -2 ~ "Малый дефицит (0 − 2% ВВП)",
      TRUE ~ "Большой дефицит (> 2% ВВП)"
    )
  ) %>%
  group_by(Режим) %>%
  summarise(
    N              = n(),
    `Средняя ИПЦ`  = round(mean(cpi,        na.rm = TRUE), 1),
    `Медиана ИПЦ`  = round(median(cpi,      na.rm = TRUE), 1),
    `Ср. ключ.ст.` = round(mean(key_rate,   na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  print(., row.names = FALSE)

cat("\n", strrep("=", 70), "\n")
cat("GMM-оценка с инструментами (бюджет как экзогенная переменная)\n")
cat(strrep("=", 70), "\n\n")

df_gmm <- df_model %>%
  arrange(year) %>%
  mutate(
    inf_expect_lag1  = lag(inf_expect, 1),
    inf_expect_lag2  = lag(inf_expect, 2),
    output_gap_lag1  = lag(output_gap, 1),
    output_gap_lag2  = lag(output_gap, 2),
    cpi_lag2         = lag(cpi, 1),
    d_reer_lag1      = lag(d_reer, 1)
  ) %>%
  filter(!is.na(inf_expect_lag2) & !is.na(output_gap_lag2))

cat("Размер выборки для GMM:", nrow(df_gmm), "\n\n")

instr_base <- ~ inf_expect_lag1 + inf_expect_lag2 +
  output_gap_lag1 + output_gap_lag2 +
  d_reer_lag1 + d_sanctions + d_covid

gmm_base <- gmm(
  g      = cpi ~ inf_expect + cpi_lag1 + output_gap,
  x      = instr_base,
  data   = df_gmm, vcov = "HAC", kernel = "Bartlett", bw = bwNeweyWest
)

gmm_fiscal <- gmm(
  g      = cpi ~ inf_expect + cpi_lag1 + output_gap + budget_def,
  x      = instr_base,        
  data   = df_gmm, vcov = "HAC", kernel = "Bartlett", bw = bwNeweyWest
)

cat("GMM базовая\n")
print(summary(gmm_base))

cat("\nGMM расширенная (+ budget_def)\n")
print(summary(gmm_fiscal))

cat("\nСравнение коэффициентов OLS vs GMM (базовая vs фискальная):\n")
coef_comp <- data.frame(
  Переменная   = names(coef(ols_base)),
  OLS_base     = round(coef(ols_base), 4),
  GMM_base     = round(coef(gmm_base)[names(coef(ols_base))], 4),
  GMM_fiscal   = round(coef(gmm_fiscal)[names(coef(ols_base))], 4)
)
print(coef_comp, row.names = FALSE)

cat(sprintf("\nGMM fiscal: budget_def коэф. = %.4f\n",
            coef(gmm_fiscal)["budget_def"]))

cat("\nПроверка силы инструментов (F первой стадии)\n")

first_stage_F <- function(endog_var, instr_formula, data) {
  fs_formula  <- reformulate(all.vars(instr_formula), response = endog_var)
  fs_model    <- lm(fs_formula, data = data)
  instr_names <- all.vars(instr_formula)
  instr_ok    <- intersect(instr_names, names(coef(fs_model)))
  ftest <- tryCatch(linearHypothesis(fs_model, instr_ok), error = function(e) NULL)
  F_stat <- if (!is.null(ftest)) round(ftest$F[2], 3) else NA
  p_val  <- if (!is.null(ftest)) round(ftest$`Pr(>F)`[2], 4) else NA
  verdict <- ifelse(!is.na(F_stat) & F_stat > 10, "сильные", "слабые")
  cat(sprintf("  %-15s | F = %7.3f | p = %.4f | %s\n",
              endog_var, ifelse(is.na(F_stat), 0, F_stat),
              ifelse(is.na(p_val), 1, p_val), verdict))
  data.frame(Эндогенная = endog_var, F_stat = F_stat, p_value = p_val,
             R2_FS = round(summary(fs_model)$r.squared, 4), Вывод = verdict)
}

fs_table <- bind_rows(
  first_stage_F("inf_expect", instr_base, df_gmm),
  first_stage_F("output_gap", instr_base, df_gmm)
)
cat("\nСводная таблица F-статистик первой стадии:\n")
print(fs_table, row.names = FALSE)
cat("\n", strrep("=", 70), "\n")
cat("GLS AR(1) — базовая и расширенная (+ бюджет)\n")
cat(strrep("=", 70), "\n\n")

gls_base <- gls(
  cpi ~ inf_expect + cpi_lag1 + output_gap + reer + d_reer + key_rate,
  data = df_model, correlation = corAR1(form = ~year), method = "ML"
)

gls_fiscal <- gls(
  cpi ~ inf_expect + cpi_lag1 + output_gap + reer + d_reer + key_rate + budget_def,
  data = df_model, correlation = corAR1(form = ~year), method = "ML"
)

cat("GLS AR(1) базовая \n");          print(summary(gls_base))
cat("\nGLS AR(1) + budget_def\n");   print(summary(gls_fiscal))

cat("\nLikelihood Ratio Test (GLS базовая и GLS + budget_def):\n")
lr_test <- anova(gls_base, gls_fiscal)
print(lr_test)

cat("\nИнтерпретация LR-теста:\n")
lr_p <- lr_test$`p-value`[2]
if (!is.na(lr_p) && lr_p < 0.05) {
  cat("Добавление budget_def ЗНАЧИМО улучшает GLS AR(1) (p < 0.05).\n")
  cat("Фискальная переменная Минфина является значимым детерминантом инфляции.\n")
} else {
  cat("Добавление budget_def не значимо улучшает модель.\n")
  cat("Возможна мультиколлинеарность с key_rate или малая выборка.\n")
}

cat("\n", strrep("=", 70), "\n")
cat("Диагностика остатков\n")
cat(strrep("=", 70), "\n\n")

resid_gls    <- residuals(gls_base,   type = "normalized")
resid_fiscal <- residuals(gls_fiscal, type = "normalized")
resid_base   <- residuals(ols_base)
resid_gmm    <- residuals(gmm_base)

lb <- function(r, lag = 5)
  round(c(stat = Box.test(r, lag = lag, type = "Ljung-Box")$statistic,
          p    = Box.test(r, lag = lag, type = "Ljung-Box")$p.value), 4)

cat("Тест Льюнга-Бокса (лаг = 5):\n")
cat(sprintf("  OLS базовая:        chi-square = %.3f, p = %.4f\n", lb(resid_base)["stat"],   lb(resid_base)["p"]))
cat(sprintf("  GLS AR(1) базовая:  chi-square = %.3f, p = %.4f\n", lb(resid_gls)["stat"],    lb(resid_gls)["p"]))
cat(sprintf("  GLS AR(1) +бюджет:  chi-square = %.3f, p = %.4f\n", lb(resid_fiscal)["stat"], lb(resid_fiscal)["p"]))
cat(sprintf("  GMM:                chi-square = %.3f, p = %.4f\n", lb(resid_gmm)["stat"],    lb(resid_gmm)["p"]))

dw_base <- dwtest(ols_base)
cat(sprintf("\nТест Дарбина-Уотсона (OLS): DW = %.4f, p = %.4f\n",
            dw_base$statistic, dw_base$p.value))

ci <- qnorm(0.975) / sqrt(nrow(df_model))
acf_df <- bind_rows(
  data.frame(lag = acf(resid_base,   plot=F, lag.max=12)$lag[-1],
             acf = acf(resid_base,   plot=F, lag.max=12)$acf[-1], model = "OLS базовая"),
  data.frame(lag = acf(resid_gls,    plot=F, lag.max=12)$lag[-1],
             acf = acf(resid_gls,    plot=F, lag.max=12)$acf[-1], model = "GLS AR(1)"),
  data.frame(lag = acf(resid_fiscal, plot=F, lag.max=12)$lag[-1],
             acf = acf(resid_fiscal, plot=F, lag.max=12)$acf[-1], model = "GLS + бюджет")
)

p_acf <- ggplot(acf_df, aes(x = lag, y = acf, fill = model)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_hline(yintercept = c(ci, -ci), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("OLS базовая" = "#4575b4",
                               "GLS AR(1)"   = "#1a9850",
                               "GLS + бюджет"= "#d73027")) +
  labs(title = "АКФ остатков (три спецификации)",
       subtitle = "Пунктиром — 95% доверительный интервал",
       x = "Лаг (лет)", y = "АКФ", fill = NULL) + theme_kp
print(p_acf)
ggsave("acf_residuals.png", p_acf, width = 10, height = 5, dpi = 150, bg = "white")

cat("\n", strrep("=", 70), "\n")
cat("Проверка робастности\n")
cat(strrep("=", 70), "\n\n")

robustness_cases <- list(
  "Полная выборка"        = df_model,
  "Без 2022"              = filter(df_model, year != 2022),
  "Без 2020"              = filter(df_model, year != 2020),
  "Без 2020 и 2022"       = filter(df_model, !year %in% c(2020, 2022)),
  "Без 1992–1993 (гипер)" = filter(df_model, year > 1993)
)

rob_tbl <- lapply(names(robustness_cases), function(nm) {
  d  <- robustness_cases[[nm]]
  m  <- lm(cpi ~ inf_expect + cpi_lag1 + output_gap + budget_def, data = d)
  hc <- coeftest(m, vcov = NeweyWest(m, lag = 2, prewhite = FALSE))
  data.frame(
    Выборка     = nm, N = nobs(m),
    `γ_f`       = round(hc["inf_expect", "Estimate"], 4),
    `γ_b`       = round(hc["cpi_lag1",   "Estimate"], 4),
    `λ`         = round(hc["output_gap", "Estimate"], 4),
    `β_budget`  = round(hc["budget_def", "Estimate"], 4),
    `p_budget`  = round(hc["budget_def", "Pr(>|t|)"], 4),
    adj_R2      = round(summary(m)$adj.r.squared, 4),
    check.names = FALSE
  )
}) %>% bind_rows()

cat("Таблица робастности (HAC, модель с бюджетным балансом):\n")
print(rob_tbl, row.names = FALSE)

cat("\nВывод по фискальному коэффициенту (β_budget) при смене выборки:\n")
rng <- range(rob_tbl$`β_budget`, na.rm = TRUE)
cat(sprintf("  Диапазон β_budget: [%.4f, %.4f]\n", rng[1], rng[2]))
if (abs(diff(rng)) < 1) {
  cat("Коэффициент стабилен при смене выборки → фискальный эффект устойчив.\n")
} else {
  cat("Значительная вариация β_budget - чувствительность к отдельным годам.\n")
}

cat("\n", strrep("=", 70), "\n")
cat("Построение графиков\n")
cat(strrep("=", 70), "\n\n")

df_plot <- df_model %>%
  mutate(
    fitted_base   = fitted(ols_base),
    fitted_gls    = fitted(gls_base),
    fitted_fiscal = fitted(models[["M6: + budget_def"]])
  )

p1 <- ggplot(df_plot, aes(x = year)) +
  geom_line(aes(y = cpi,           color = "Факт"),               linewidth = 1.3) +
  geom_line(aes(y = fitted_base,   color = "Гибридная PC (OLS)"), linewidth = 1.0, linetype = "dashed") +
  geom_line(aes(y = fitted_fiscal, color = "PC + бюджет (OLS)"),  linewidth = 1.0, linetype = "dotdash") +
  geom_line(aes(y = fitted_gls,    color = "GLS AR(1)"),          linewidth = 1.0, linetype = "dotted") +
  scale_color_manual(values = c(
    "Факт"                  = "#d73027",
    "Гибридная PC (OLS)"    = "#4575b4",
    "PC + бюджет (OLS)"     = "#f46d43",
    "GLS AR(1)"             = "#1a9850"
  )) +
  labs(title    = "Инфляция ИПЦ: факт и оценки моделей кривой Филлипса",
       subtitle = "Россия, 2006–2026 гг. (включая фискальную спецификацию)",
       x = "Год", y = "Инфляция ИПЦ (%)", color = NULL) + theme_kp
print(p1)
ggsave("fig1_fact_vs_models.png", p1, width = 11, height = 5, dpi = 150, bg = "white")

resid_df <- df_plot %>%
  mutate(
    resid_base   = cpi - fitted_base,
    resid_fiscal = cpi - fitted_fiscal
  ) %>%
  select(year, resid_base, resid_fiscal) %>%
  pivot_longer(-year, names_to = "model", values_to = "resid") %>%
  mutate(model = dplyr::recode(model,
                               resid_base   = "OLS базовая (без бюджета)",
                               resid_fiscal = "OLS + budget_def (Минфин)"
  ))

p2 <- ggplot(resid_df, aes(x = year, y = resid, fill = resid > 0)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.7) +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#1a9850", "FALSE" = "#d73027")) +
  labs(title    = "Остатки: базовая модель и модель с бюджетным балансом",
       subtitle = "Зелёный — переоценка модели, красный — недооценка",
       x = "Год", y = "Остаток (п.п.)") + theme_kp
print(p2)
ggsave("fig2_residuals.png", p2, width = 11, height = 6, dpi = 150, bg = "white")

p3a <- ggplot(df_model, aes(x = year)) +
  geom_col(aes(y = budget_def, fill = budget_def > 0), alpha = 0.75, width = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.7, color = "grey30") +
  scale_fill_manual(values = c("TRUE" = "#1a9850", "FALSE" = "#d73027"),
                    labels = c("Дефицит", "Профицит"),
                    name   = NULL) +
  labs(title = "Бюджетный баланс Минфина РФ (% ВВП)",
       x = "Год", y = "% ВВП") + theme_kp

p3b <- ggplot(df_model, aes(x = year)) +
  geom_line(aes(y = cpi),       color = "#d73027", linewidth = 1.2) +
  geom_line(aes(y = key_rate),  color = "#4575b4", linewidth = 1.0, linetype = "dashed") +
  labs(title = "Инфляция ИПЦ (красный) и ключевая ставка ЦБ (синий)",
       x = "Год", y = "%") + theme_kp

p3 <- p3a / p3b +
  plot_annotation(
    title    = "Фискальная политика Минфина, монетарная политика ЦБ и инфляция",
    subtitle = "Россия, 2006–2026 гг.",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )
print(p3)
ggsave("fig3_fiscal_monetary.png", p3, width = 11, height = 7, dpi = 150, bg = "white")

df_scatter <- df_model %>%
  mutate(Эпоха = case_when(
    year <= 2007 ~ "2000–2007\n(нефтяной бум)",
    year <= 2013 ~ "2008–2013\n(кризис + восст.)",
    year <= 2021 ~ "2014–2021\n(санкции/ковид)",
    TRUE         ~ "2022–2026\n(Санкции)"
  ))

p4 <- ggplot(df_scatter, aes(x = budget_def, y = cpi, color = Эпоха, label = year)) +
  geom_point(size = 3.5) +
  geom_text_repel(size = 2.7, max.overlaps = 20) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  scale_color_manual(values = c("#d73027","#f46d43","#4575b4","#1a9850","#756bb1")) +
  labs(
    title    = "Бюджетный баланс Минфина РФ и инфляция ИПЦ",
    subtitle = "Каждая точка — один год; пунктир — общий линейный тренд",
    x = "Бюджетный баланс (% ВВП, + = профицит, − = дефицит)",
    y = "Инфляция ИПЦ (%)", color = "Эпоха"
  ) + theme_kp
print(p4)
ggsave("fig4_budget_vs_inflation.png", p4, width = 10, height = 7, dpi = 150, bg = "white")

rob_long <- rob_tbl %>%
  select(Выборка, `γ_f`, `γ_b`, `λ`, `β_budget`) %>%
  pivot_longer(-Выборка, names_to = "Коэф.", values_to = "Оценка")

p5 <- ggplot(rob_long, aes(x = Выборка, y = Оценка, fill = `Коэф.`)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("γ_f"="royalblue","γ_b"="#d73027",
                               "λ"="#1a9850","β_budget"="#f46d43")) +
  labs(title    = "Стабильность коэффициентов при исключении отдельных лет",
       subtitle = "Модель: cpi ~ inf_expect + cpi_lag1 + output_gap + budget_def (HAC)",
       x = NULL, y = "Оценка коэффициента", fill = "Коэффициент") +
  theme_kp + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 10))
print(p5)
ggsave("fig5_robustness_budget.png", p5, width = 11, height = 5, dpi = 150, bg = "white")

cat("\n", strrep("=", 70), "\n")
cat("ИТОГОВОЕ СРАВНЕНИЕ МОДЕЛЕЙ\n")
cat(strrep("=", 70), "\n\n")

final_tbl <- data.frame(
  Модель        = c("OLS базовая (HAC)",
                    "OLS + budget_def (HAC)",
                    "GLS AR(1) базовая",
                    "GLS AR(1) + budget_def",
                    "GMM базовая",
                    "GMM + budget_def"),
  N             = c(nobs(ols_base),
                    nobs(models[["M6: + budget_def"]]),
                    nobs(gls_base),
                    nobs(gls_fiscal),
                    nrow(df_gmm), nrow(df_gmm)),
  R2_adj        = c(round(summary(ols_base)$adj.r.squared, 4),
                    round(summary(models[["M6: + budget_def"]])$adj.r.squared, 4),
                    NA, NA, NA, NA),
  AIC           = c(round(AIC(ols_base), 2),
                    round(AIC(models[["M6: + budget_def"]]), 2),
                    round(AIC(gls_base), 2),
                    round(AIC(gls_fiscal), 2),
                    NA, NA),
  LjBox_p5      = c(round(lb(resid_base)["p"],   4),
                    NA,
                    round(lb(resid_gls)["p"],    4),
                    round(lb(resid_fiscal)["p"], 4),
                    round(lb(resid_gmm)["p"],    4),
                    NA),
  beta_budget   = c(NA,
                    round(hac_fiscal["budget_def","Estimate"], 4),
                    NA,
                    round(summary(gls_fiscal)$tTable["budget_def","Value"], 4),
                    NA,
                    round(coef(gmm_fiscal)["budget_def"], 4)),
  p_budget      = c(NA,
                    round(hac_fiscal["budget_def","Pr(>|t|)"], 4),
                    NA,
                    round(summary(gls_fiscal)$tTable["budget_def","p-value"], 4),
                    NA,
                    round(summary(gmm_fiscal)$coef["budget_def","Pr(>|t|)"], 4))
)
print(final_tbl, row.names = FALSE)

cat("\nСводка: стационарность\n")
print(ur_table[ , c("Переменная","ADF_вывод","KPSS_вывод")], row.names = FALSE)

cat("\nСила инструментов GMM\n")
print(fs_table[ , c("Эндогенная","F_stat","Вывод")], row.names = FALSE)

cat("\nСохранённые графики: \n",
    "  fig0_stationarity_dynamics.png — динамика рядов (вкл. бюджет)\n",
    "  fig1_fact_vs_models.png        — факт. и модели (базовая + фискальная)\n",
    "  fig2_residuals.png             — остатки двух спецификаций\n",
    "  fig3_fiscal_monetary.png       — бюджет Минфина + ставка ЦБ + инфляция\n",
    "  fig4_budget_vs_inflation.png   — scatter: баланс vs инфляция по эпохам\n",
    "  fig5_robustness_budget.png     — стабильность β_budget\n",
    "  acf_residuals.png              — АКФ трёх моделей\n")