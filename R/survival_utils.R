plot_kaplan_meier <- function(event, time_to_event, group, colors = NULL,
                               time_range = c(0, 120), title = NULL,
                               ylab = "Recurrence free survival",
                               censor_shape = "+") {
  library(survival)
  library(survminer)

  time_to_event[time_to_event == 0] <- NA
  df <- data.frame(
    event = as.numeric(event),
    time_to_event = as.numeric(time_to_event),
    group = group
  )
  df <- df[complete.cases(df), ]

  fit <- survfit(Surv(time_to_event, event) ~ group, data = df)
  names(fit$strata) <- gsub("group=", "", names(fit$strata))
  names(colors) <- NULL

  ggsurvplot(
    fit, data = df,
    pval = TRUE,
    risk.table = TRUE, risk.table.col = "strata",
    surv.median.line = "none",
    ggtheme = theme_test(),
    palette = colors,
    xlim = time_range,
    break.time.by = 24,
    tables.y.text = TRUE,
    censor.shape = censor_shape,
    legend = "top",
    title = title,
    ylab = ylab
  )
}

#' Truncate a recurrence-free-survival cohort at `n` days: events occurring
#' after the horizon are recoded as censored at the horizon 
truncate_survival <- function(event, time_to_event, horizon_days) {
  time_to_event <- as.numeric(time_to_event)
  event <- as.numeric(event)
  past_horizon <- time_to_event > horizon_days
  event[past_horizon] <- 0
  time_to_event[past_horizon] <- horizon_days
  list(event = event, time_to_event = time_to_event)
}
