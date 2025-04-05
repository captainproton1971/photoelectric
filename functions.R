# functions.R



multimeter_file <- "./metex_M3800.yaml"
# Load the YAML configuration file
dmm_config <- yaml::yaml.load_file(multimeter_file)

voltmeter_config <- as.data.frame(dmm_config$voltmeter_V)

tidy.gls <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  est <- coef(x)
  se  <- sqrt(diag(vcov(x)))
  df  <- x$dims$N - length(est)
  res <- data.frame(term = names(est), estimate = est, std.error = se)

  if (conf.int) {
    alpha <- 1 - conf.level
    tval <- qt(1 - alpha / 2, df)
    res$conf.low <- est - tval * se
    res$conf.high <- est + tval * se
  }

  res
}

smart_format <- function(x, sci_thresh = 1e-3, digits = 4, sci_digits = 3) {
  x <- as.numeric(x)  # ensure numeric
  sapply(x, function(val) {
    if (is.na(val)) return(NA)
    if (abs(val) < sci_thresh) {
      formatC(val, format = "e", digits = sci_digits)
    } else {
      round(val, digits)
    }
  })
}

pretty_uncert <- function(x, dx){
  # A small function for displaying a value and uncertainty to appropriate number of sig. figs.  General rule, keep one uncertain digit
  # unless leading digit is 1 in which case, retain two uncertain digtis.

  pos <- -floor(log10(abs(dx)))
  u_digits <- dx * 10^pos

  if (round(u_digits < 2)){
    pos <- pos + 1
  }
  f_string <- paste0("%.", pos, "f")

  x_str <- sprintf(f_string, round(x*10^pos)/10^pos)
  u_str <- sprintf(f_string, round(dx*10^pos)/10^pos)
  return(x_str)
  #return(paste0(x_str, "(", u_str,")"))
}

# Function for Iterated GLS Regression
iterated_gls <- function(df, include_constant = TRUE, tolerance = 1e-3, max_iter = 100) {
  # Initial slope estimate (OLS)
  beta_old <- 0.1
  beta_new <- 1
  iter <- 0

  while(abs((beta_new - beta_old)/beta_old) > tolerance && iter < max_iter) {
    iter <- iter + 1
    beta_old <- beta_new

    # Compute effective y variance
    df$adjusted_vars <- df$sigma_y^2 + beta_old^2 * df$sigma_x^2

    # Fit GLS with updated weights
    if(include_constant) {
      model <- gls(y ~ x, data = df, weights = varFixed(~ adjusted_vars))
    } else {
      model <- gls(y ~ x - 1, data = df, weights = varFixed(~ adjusted_vars))
    }

    beta_new <- as.numeric(coef(model)[["x"]])
  }

  list(model = model, iterations = iter, converged = (iter < max_iter))
}


create_plot <- function(df, model, units = "eV") {
  # Rescale x-axis
  x_scale_factor <- 1e14
  df$x <- df$x / x_scale_factor
  df$sigma_x <- df$sigma_x / x_scale_factor

  # Compute y-scale factor (autoscale to bring values into 1–10 range)
  max_y <- max(abs(df$y))
  scale_exp <- floor(log10(max_y))
  y_scale_factor <- if (scale_exp < -1 || scale_exp > 2) 10^scale_exp else 1
  df$y <- df$y / y_scale_factor
  df$sigma_y <- df$sigma_y / y_scale_factor

  # x values for prediction
  x_max_extended <- 1.1 * max(df$x)
  x_vals <- seq(0.0, x_max_extended, length.out = 200)
  new_df <- data.frame(x = x_vals)

  # Model predictions (x back to Hz for prediction)
  X <- model.matrix(~ x, data = data.frame(x = x_vals * x_scale_factor))
  fit <- as.vector(X %*% coef(model))

  # Prediction errors
  vcov_beta <- vcov(model)
  se_fit <- sqrt(diag(X %*% vcov_beta %*% t(X)))
  dof <- model$dims$N - length(coef(model))
  t_crit <- qt(0.975, df = dof)
  upr <- fit + t_crit * se_fit
  lwr <- fit - t_crit * se_fit

  # Rescale predictions for plotting
  pred_df <- data.frame(x = x_vals, fit = fit / y_scale_factor,
                        lwr = lwr / y_scale_factor, upr = upr / y_scale_factor)

  # Build plot
  p <- ggplot() +
    geom_line(data = pred_df, aes(x = x, y = fit), color = "blue") +
    geom_ribbon(data = pred_df, aes(x = x, ymin = lwr, ymax = upr),
                alpha = 0.2, fill = "blue") +
    geom_point(data = df, aes(x = x, y = y)) +
    geom_errorbar(data = df,
                  aes(x = x, ymin = y - sigma_y, ymax = y + sigma_y,
                      width = sigma_x),
                  color = "black") +
    geom_errorbarh(data = df,
                   aes(y = y, xmin = x - sigma_x, xmax = x + sigma_x,
                       height = sigma_y),
                   color = "black") +
    expand_limits(x = c(0, x_max_extended)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    theme_minimal() +
    labs(
      x = expression("Frequency (" * 10^{14} * " Hz)"),
      y = if (y_scale_factor == 1) {
        paste0("Electron Kinetic Energy (", units, ")")
      } else {
        bquote("Electron Kinetic Energy (" * 10^.(scale_exp) * " " * .(units) * ")")
      }
    )

  return(p)
}








#' Get Rounded Values
#'
#' Helper function for rounding a reading and
#' associated uncertainty  measurements. It rounds the reading based on
#' the provided uncertainty and  decimals, ensuring that the rounding
#' reflects the significant figures  dictated by the uncertainty.

#' @param x A numeric vector of length 3: reading, uncertainty, and
#' 	number of    decimal places.
#'
#' @return A list of strings of length 4 corresponding to (1) the
#' 	rounded  value, (2) the rounded uncertainty, (3)
#' 	rounded_value \eqn{\pm} rounded_uncertainty,  and (4)
#' 	rounded_value(uncertain_digits)
#'
#' @examples output <-
#'   get_rounded(c(123.456, 0.789, 2))
#' print(output)
# @export
get_rounded <- function(x) {
  reading <- x[1]
  uncertainty <- x[2]
  decimals <- x[3]

  exponent <- floor(log10(abs(uncertainty)))
  mantissa <- uncertainty / 10^exponent

  round_val <- -exponent
  if (round(mantissa, 0) < 2) {
    round_val <- round_val + 1
    mantissa <- mantissa * 10
  }

  rounded_reading <- format(
    round(reading, min(round_val, decimals)),
    nsmall = min(round_val, decimals)
  )
  rounded_uncertainty <- format(round(uncertainty, min(round_val, decimals)),
                                nsmall = min(round_val, decimals)
  )

  return(
    c(
      rounded_reading,
      rounded_uncertainty,
      paste0(rounded_reading, " \u00B1 ", rounded_uncertainty),
      paste0(rounded_reading, "(", round(mantissa,0), ")")
    )
  )
}

#' Add Errors in Quadrature
#
#' This function computes the combined standard uncertainty for a set
#' of measurements  with individual uncertainties, assuming that the
#' errors are uncorrelated. It uses  the method of adding errors in
#' quadrature, which is a standard procedure in error  analysis for
#' combining multiple independent error terms.
#
#' @param x A numeric vector of uncertainties for each measurement.
#' @return Returns the combined standard uncertainty as a numeric
#' 	value.
#' @examples
#' error_values <- c(0.1, 0.2, 0.05)
#' combined_error <- error_add(error_values)
#' @export
error_add <- function(x) {
  # Simple adding by quadratures for error propagation
  # x is an array of numerics.
  # e.g., x=c(3,4,5,6,7)
  # (e.g., continued) value <- error_add(x)
  # value =
  if (!is.numeric(x)) {
    return("Not numeric")
  } else {
    return(sqrt(sum(x * x)))
  }
}

#' count_decimal_places
#' This helper function counts the decimal places of a floating-point
#' number passed as a sting `input_string.` This format is required so
#' as to not lose decimal places when passing a value such as 1.730.
#
#' @param input_string A string representation of a floating point
#' 	number.
#' @return Returns the number of decimal places in the input_string.
#' @examples
#' x <- "1.302"
#' x_decimals <- count_decimal_places(x)
#' print(x_decimals)
count_decimal_places <- function(input_string) {
  # Check if the string contains a decimal point.
  if (grepl("\\.", input_string)) {
    # Split the string at the decimal point.
    parts <- unlist(strsplit(input_string, "\\.", fixed = FALSE))
    # The number of decimal places is the length of the substring after
    # the decimal point.
    num_decimals <- nchar(parts[2])
  } else {
    # If theres no decimal point, the number of decimal places is 0.
    num_decimals <- 0
  }

  return(num_decimals)
}
#' Determine Suitable Index for Measurement
#'
#' Determines the suitable index based on the measurement scale and resolution.
#' For 'auto' scale, it finds indices based on range and resolution conditions.
#' For specific scales, it matches the scale name directly.
#'
#' @param scale The measurement scale to use, or 'auto' to determine based on
#'   conditions.
#' @param multiplier The multiplier to adjust measurement scale.
#' @param config_df Data frame with meter configuration details.
#' @param x_val The measurement value.
#' @param emp_rez The empirical resolution calculated from the measurement.
#' @return Integer index indicating the suitable configuration row, or NULL if
#'   not found.
determine_idx <- function(scale, multiplier, config_df, x_val, emp_rez) {
  if (scale == "auto") {
    idx1 <- which(
      config_df$range * config_df$unit_multiplier >= x_val * multiplier
    )
    idx2 <- which(config_df$resolution >= emp_rez)
    if (length(idx1) == 0 || length(idx2) == 0) {
      return(NULL)
    } # Indicate no suitable index found
    return(max(idx1[1], idx2[1]))
  } else {
    idx <- which(config_df$scale_names == scale)
    if (length(idx) == 0) {
      return(NULL)
    } # Indicate no suitable index found
    return(idx[1])
  }
}



#' Multimeter Uncertainty Calculations
#'
#' Description of how `meter_uncertainty` works and its role as a core function for calculating
#' measurement uncertainty for various types of meters.
#'
#' @name meter_uncertainty
#' @param input_string A string containing the meter reading, specific to the type of meter.
#' @param multiplier The measurement scale multiplier.
#' @param scale The measurement scale, or 'auto' for automatic determination.
#' @param config_df A data frame containing the meter's configuration and accuracy specs.
#' @return A list with calculated measurement uncertainty details.
meter_uncertainty <- function(input_string, multiplier, scale, config_df) {
  original_options <- options("warn", "show.error.messages")
  on.exit(options(original_options))

  tryCatch(
    {
      suppressWarnings({
        x_val <- as.numeric(input_string)
        if (is.na(x_val)) {
          result <- rep(paste0(input_string, " is not castable to float."), 4)
          return(result)
        }

        decimals <- count_decimal_places(input_string)
        emp_rez <- 10^-(decimals - log10(multiplier))

        # Simplified idx determination
        idx <- determine_idx(scale, multiplier, config_df, x_val, emp_rez)
        if (is.null(idx)) {
          return(rep(paste0("Improper scale `scale = ", scale, "."), 4))
        }

        res_decimals <- -log10(abs(config_df$resolution[idx] / multiplier))

        # Find the value of uncertainty
        uncertainty <- config_df$multiplier[idx] * x_val +
          config_df$last_digit[idx] * 10^(-res_decimals)

        result <- get_rounded(c(x_val, uncertainty, decimals))
      })
    },
    error = function(e) {
      return(
        rep(paste0("An error occurred during calculation: ", e$message), 4)
      )
    }
  )

  return(result)
}

voltmeter_uncertainty <- function(input_string, unit = "V", scale = "auto") {
  multiplier <- switch(unit,
                       "mV" = 1.0e-03,
                       "V" = 1.0
  )
  if (is.null(multiplier)) {
    result <- rep(paste0("Invalid unit ", unit, "."), 4)
    return(result)
  }

  config_df <- voltmeter_config

  result <- meter_uncertainty(input_string, multiplier, scale, config_df)
  return(result)
}

process_voltages <- function(df_numeric, element_function = voltmeter_uncertainty) {
  # Step 1: Apply the element-wise function to each stringified entry
  results <- lapply(df_numeric, function(col) {
    lapply(col, function(x) element_function(as.character(x)))
  })

  # Helper to extract nth element from each result
  extract_nth_output <- function(n) {
    suppressWarnings({
    extracted <- lapply(results, function(col) sapply(col, `[`, n))
    as.data.frame(lapply(extracted, as.numeric))
    })
  }

  # Step 2: Extract values and uncertainties
  values <- extract_nth_output(1)
  uncertainties <- extract_nth_output(2)

  # Step 3: Compute stats

  row_counts <- apply(df_numeric, 1, function(x) sum(!is.na(x)))
  means <- rowMeans(values, na.rm = TRUE)

  row_var_vals <- apply(values, 1, function(v) {
    if (sum(!is.na(v)) <= 1) 0 else var(v, na.rm = TRUE)
  })

  uncert_var_vals <- apply(uncertainties^2, 1, function(v) {
    mean(v, na.rm = TRUE)
  })

  row_vars <- row_var_vals + uncert_var_vals
  sigmas <- sqrt(row_vars / row_counts)
  sigmas[row_counts == 0] <- NA  # optional safety check

  return(data.frame(y = means, sigma_y = sigmas))
}

