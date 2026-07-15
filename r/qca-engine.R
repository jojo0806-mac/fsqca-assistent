suppressPackageStartupMessages({
  library(QCA)
  library(jsonlite)
})

read_qca_input <- function() {
  fromJSON("/qca-input.json", simplifyVector = TRUE, simplifyDataFrame = TRUE)
}

write_qca_output <- function(value) {
  toJSON(value, auto_unbox = TRUE, dataframe = "rows", na = "null", null = "null", digits = 15)
}

direct_calibrate <- function(x, anchor) {
  full <- as.numeric(anchor$full)
  crossover <- as.numeric(anchor$crossover)
  non <- as.numeric(anchor$non)
  if (identical(anchor$direction, "negative")) {
    x <- -x
    old_full <- full
    full <- -non
    non <- -old_full
    crossover <- -crossover
  }
  if (!(full > crossover && crossover > non)) {
    stop("校准锚点必须满足：完全隶属 > 交叉点 > 完全不隶属")
  }
  logit95 <- log(0.95 / 0.05)
  scaled <- ifelse(
    x >= crossover,
    ((x - crossover) / (full - crossover)) * logit95,
    ((x - crossover) / (crossover - non)) * logit95
  )
  membership <- 1 / (1 + exp(-scaled))
  membership[membership == 0.5] <- 0.500001
  pmin(0.999, pmax(0.001, membership))
}

qca_calibrate <- function() {
  input <- read_qca_input()
  dat <- input$rows
  variables <- as.character(input$variables)
  missing_mode <- input$missing
  numeric_data <- lapply(variables, function(variable) {
    suppressWarnings(as.numeric(as.character(dat[[variable]])))
  })
  names(numeric_data) <- variables
  missing_rows <- Reduce(`|`, lapply(numeric_data, is.na))
  removed <- if (identical(missing_mode, "drop")) sum(missing_rows) else 0
  imputed <- 0

  if (identical(missing_mode, "drop")) {
    keep <- which(!missing_rows)
  } else {
    keep <- seq_len(nrow(dat))
    for (variable in variables) {
      values <- numeric_data[[variable]]
      replacement <- if (identical(missing_mode, "median")) {
        median(values, na.rm = TRUE)
      } else {
        mean(values, na.rm = TRUE)
      }
      imputed <- imputed + sum(is.na(values))
      values[is.na(values)] <- replacement
      numeric_data[[variable]] <- values
    }
  }

  fuzzy <- lapply(variables, function(variable) {
    values <- numeric_data[[variable]][keep]
    if (identical(input$modes[[variable]], "keep")) {
      if (any(values < 0 | values > 1)) {
        stop(sprintf("变量‘%s’选择了保留原值，但存在超出 0–1 范围的数值", variable))
      }
      values
    } else {
      direct_calibrate(values, input$anchors[[variable]])
    }
  })
  names(fuzzy) <- variables

  rows <- lapply(seq_along(keep), function(i) {
    list(
      rowIndex = keep[[i]] - 1,
      caseId = as.character(dat[[input$caseColumn]][keep[[i]]]),
      fuzzy = as.list(vapply(fuzzy, function(column) column[[i]], numeric(1)))
    )
  })
  write_qca_output(list(rows = rows, removed = removed, imputed = imputed))
}

analysis_data <- function(input) {
  dat <- input$data
  variables <- c(as.character(input$conditions), input$outcome)
  for (variable in variables) dat[[variable]] <- as.numeric(dat[[variable]])
  dat
}

qca_necessity <- function() {
  input <- read_qca_input()
  dat <- analysis_data(input)
  outcome_values <- dat[[input$outcome]]
  result <- list()
  index <- 1
  for (condition in as.character(input$conditions)) {
    for (negated in c(FALSE, TRUE)) {
      condition_values <- if (negated) 1 - dat[[condition]] else dat[[condition]]
      intersection <- sum(pmin(condition_values, outcome_values))
      result[[index]] <- list(
        condition = condition,
        negated = negated,
        consistency = intersection / sum(outcome_values),
        coverage = intersection / sum(condition_values)
      )
      index <- index + 1
    }
  }
  write_qca_output(result)
}

make_truth_table <- function(input) {
  dat <- analysis_data(input)
  labels <- as.character(input$conditions)
  safe_conditions <- paste0("C", seq_along(labels))
  qca_data <- data.frame(lapply(labels, function(label) dat[[label]]), check.names = FALSE)
  names(qca_data) <- safe_conditions
  qca_data$OUTCOME <- dat[[input$outcome]]
  qca_data$CASE_ID <- dat$`_case_id`
  tt <- truthTable(
    qca_data,
    outcome = "OUTCOME",
    conditions = safe_conditions,
    incl.cut = as.numeric(input$thresholds$consistency),
    n.cut = as.numeric(input$thresholds$frequency),
    pri.cut = as.numeric(input$thresholds$pri),
    complete = TRUE,
    show.cases = TRUE
  )
  list(tt = tt, safe_conditions = safe_conditions, labels = labels)
}

qca_truth_table <- function() {
  input <- read_qca_input()
  built <- make_truth_table(input)
  tt <- built$tt
  table <- tt$tt
  conditions <- built$safe_conditions
  present <- which(table$n > 0)
  table <- table[present, , drop = FALSE]
  out_rank <- ifelse(table$OUT == "1", 0, ifelse(table$OUT == "0", 1, 2))
  pri_value <- suppressWarnings(as.numeric(table$PRI))
  incl_value <- suppressWarnings(as.numeric(table$incl))
  table <- table[order(out_rank, -pri_value, -incl_value, -table$n), , drop = FALSE]

  rows <- lapply(seq_len(nrow(table)), function(i) {
    bits <- as.integer(table[i, conditions, drop = TRUE])
    code <- sum(bits * (2 ^ rev(seq_along(bits) - 1)))
    list(
      bits = bits,
      code = code,
      frequency = as.integer(table$n[[i]]),
      consistency = as.numeric(table$incl[[i]]),
      pri = as.numeric(table$PRI[[i]]),
      outcome = if (table$OUT[[i]] == "1") 1 else 0,
      cases = if (nzchar(table$cases[[i]])) strsplit(table$cases[[i]], ",", fixed = TRUE)[[1]] else character(),
      remainder = table$OUT[[i]] == "?"
    )
  })
  write_qca_output(rows)
}

term_pattern <- function(term, conditions) {
  literals <- strsplit(term, "*", fixed = TRUE)[[1]]
  paste(vapply(conditions, function(condition) {
    if (condition %in% literals) "1" else if (paste0("~", condition) %in% literals) "0" else "-"
  }, character(1)), collapse = "")
}

display_term <- function(term, conditions, labels) {
  literals <- strsplit(term, "*", fixed = TRUE)[[1]]
  paste(vapply(literals, function(literal) {
    negated <- startsWith(literal, "~")
    bare <- if (negated) substring(literal, 2) else literal
    index <- match(bare, conditions)
    label <- if (is.na(index)) bare else labels[[index]]
    paste0(if (negated) "~" else "", label)
  }, character(1)), collapse = "*")
}

normal_solution <- function(solution, kind, conditions, labels, intermediate = FALSE) {
  if (intermediate) {
    selected <- solution$i.sol[[1]]
    terms <- selected$solution
    ic <- selected$IC
  } else {
    terms <- solution$solution[[1]]
    ic <- solution$IC
  }
  fit <- if (!is.null(ic$overall)) {
    ic$overall$sol.incl.cov[1, , drop = FALSE]
  } else if (!is.null(ic$individual) && length(ic$individual) > 0) {
    ic$individual[[1]]$sol.incl.cov[1, , drop = FALSE]
  } else {
    ic$sol.incl.cov[1, , drop = FALSE]
  }
  path_fit <- if (!is.null(ic$individual) && length(ic$individual) > 0) {
    ic$individual[[1]]$incl.cov
  } else {
    ic$incl.cov
  }
  display_terms <- vapply(terms, display_term, character(1), conditions = conditions, labels = labels)
  patterns <- vapply(terms, term_pattern, character(1), conditions = conditions)
  metric_value <- function(row, column) {
    if (is.null(path_fit) || nrow(path_fit) < row || !column %in% colnames(path_fit)) return(NA_real_)
    as.numeric(path_fit[[column]][[row]])
  }
  path_metrics <- lapply(seq_along(terms), function(index) {
    list(
      term = unname(display_terms[[index]]),
      pattern = unname(patterns[[index]]),
      consistency = metric_value(index, "inclS"),
      pri = metric_value(index, "PRI"),
      rawCoverage = metric_value(index, "covS"),
      uniqueCoverage = metric_value(index, "covU")
    )
  })
  list(
    kind = kind,
    patterns = patterns,
    terms = unname(display_terms),
    formula = paste(display_terms, collapse = " + "),
    consistency = as.numeric(fit$inclS[[1]]),
    pri = as.numeric(fit$PRI[[1]]),
    coverage = as.numeric(fit$covS[[1]]),
    pathMetrics = unname(path_metrics)
  )
}

qca_solutions <- function() {
  input <- read_qca_input()
  built <- make_truth_table(input)
  tt <- built$tt
  if (!any(tt$tt$OUT == "1")) stop("当前阈值下没有结果为 1 的组态")
  conditions <- built$safe_conditions
  labels <- built$labels
  complex <- minimize(tt, include = "", details = TRUE, show.cases = TRUE)
  parsimonious <- minimize(tt, include = "?", details = TRUE, show.cases = TRUE)
  directions <- vapply(seq_along(conditions), function(index) {
    direction <- input$directions[[labels[[index]]]]
    if (identical(direction, "positive")) "1" else if (identical(direction, "negative")) "0" else "-"
  }, character(1))
  intermediate <- minimize(
    tt,
    include = "?",
    dir.exp = paste(directions, collapse = ","),
    details = TRUE,
    show.cases = TRUE
  )
  intermediate_result <- if (all(directions == "-")) {
    normal_solution(complex, "intermediate", conditions, labels)
  } else {
    normal_solution(intermediate, "intermediate", conditions, labels, intermediate = TRUE)
  }
  write_qca_output(list(
    normal_solution(complex, "complex", conditions, labels),
    intermediate_result,
    normal_solution(parsimonious, "parsimonious", conditions, labels)
  ))
}
