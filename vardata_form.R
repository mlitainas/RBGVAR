vardata_form = function(data , lags = 1, const = TRUE, trend = FALSE, trend_qua = FALSE, star = NULL, 
                        xlags = 0, det_controls = NULL, ex= NULL, ex_lag = NULL,  feedback_only = FALSE ){
  
  # -------------------------------------------------------------------------- #
  # Helper functions
  # -------------------------------------------------------------------------- #
  check_logical_scalar <- function(x, argument_name) {
    if (
      length(x) != 1L ||
      !is.logical(x) ||
      is.na(x)
    ) {
      stop(
        sprintf("'%s' must be either TRUE or FALSE.", argument_name),
        call. = FALSE
      )
    }
  }
  
  check_integer_scalar <- function(x, argument_name, minimum = 0L) {
    if (
      length(x) != 1L ||
      !is.numeric(x) ||
      is.na(x) ||
      !is.finite(x) ||
      x < minimum ||
      x %% 1 != 0
    ) {
      description <- if (minimum == 0L) {
        "a non-negative integer"
      } else if (minimum == 1L) {
        "a positive integer"
      } else {
        paste0("an integer greater than or equal to ", minimum)
      }
      stop(sprintf("'%s' must be %s.", argument_name, description), call. = FALSE)
    }
  }
  
  check_numeric_block <- function(x, argument_name, expected_rows) {
    if (!(is.data.frame(x) || is.matrix(x))) {
      stop(
        sprintf("'%s' must be a data frame, tibble, or matrix.", argument_name),
        call. = FALSE
      )
    }
    
    if (ncol(x) < 1L) {
      stop(
        sprintf("'%s' must contain at least one column.", argument_name),
        call. = FALSE
      )
    }
    
    if (nrow(x) != expected_rows) {
      stop(
        sprintf(
          "'%s' must have the same number of rows as 'data'.",
          argument_name
        ),
        call. = FALSE
      )
    }
    
    if (!all(vapply(as.data.frame(x), is.numeric, logical(1)))) {
      stop(
        sprintf("All columns of '%s' must be numeric.", argument_name),
        call. = FALSE
      )
    }
    
    if (anyNA(x)) {
      stop(
        sprintf("'%s' contains missing values.", argument_name),
        call. = FALSE
      )
    }
    
    if (!all(is.finite(as.matrix(x)))) {
      stop(
        sprintf("'%s' contains non-finite values.", argument_name),
        call. = FALSE
      )
    }
    
    if (is.null(colnames(x)) || any(colnames(x) == "")) {
      stop(
        sprintf("All columns of '%s' must have names.", argument_name),
        call. = FALSE
      )
    }
    
    if (anyDuplicated(colnames(x))) {
      stop(
        sprintf("Column names in '%s' must be unique.", argument_name),
        call. = FALSE
      )
    }
  }
  
  
  # -------------------------------------------------------------------------- #
  # Check the consistency of the inputs
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # 1. Check data
  # -------------------------------------------------------------------------- #
  
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or tibble.",call. = FALSE)
  }
  
  if (nrow(data) < 1L) {
    stop("'data' must contain at least one observation.", call. = FALSE )
  }
  
  if (ncol(data) < 1L) {
    stop("'data' must contain at least one variable.",call. = FALSE)
  }
  
  if (
    is.null(names(data)) ||
    anyNA(names(data)) ||
    any(names(data) == "")
  ) {
    stop(
      "All columns of 'data' must have names.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(names(data))) {
    stop("Column names in 'data' must be unique.",call. = FALSE)
  }
  
  # Identify and remove the time-identifier column.
  time_id_name <- names(data)[1L]
  
  if (
    time_id_name %in% c("Date", "date", "Time", "time")
  ) {
    
    if (anyNA(data[[1L]])) {
      stop(
        "The time identifier contains missing values.",
        call. = FALSE
      )
    }
    
    if (anyDuplicated(data[[1L]])) {
      stop(
        "The time identifier contains duplicate values.",
        call. = FALSE
      )
    }
    
    TimeID <- data.frame(TimeID = data[[1L]])
    
    data <- data[-1L]
    
    if (ncol(data) < 1L) {
      stop(
        paste0(
          "'data' must contain at least one endogenous ",
          "variable in addition to the time identifier."
        ),
        call. = FALSE
      )
    }
    
    message("Time identifier detected: ", time_id_name )
    
  } else {
    
    TimeID <- NULL
    message("No time identifier detected.")
    
  }
  
  # Check the endogenous variables after removing the time identifier.
  if (!all(vapply(data, is.numeric, logical(1)))) {
    stop(
      "All endogenous variables in 'data' must be numeric.",
      call. = FALSE
    )
  }
  
  if (anyNA(data)) {
    stop("'data' contains missing values.",call. = FALSE)
  }
  
  if (!all(is.finite(as.matrix(data)))) {
    stop("'data' contains non-finite values.",call. = FALSE)
  }
  
  # Number of observations after removing the time identifier.
  number_of_observations <- nrow(data)

  # -------------------------------------------------------------------------- #
  # 2. Check lags
  # -------------------------------------------------------------------------- #
  check_integer_scalar(lags,"lags", minimum = 1L)
  
  # -------------------------------------------------------------------------- #
  # 3. Check const
  # -------------------------------------------------------------------------- #
  check_logical_scalar(const,"const")
  
  # -------------------------------------------------------------------------- #
  # 4. Check trend
  # -------------------------------------------------------------------------- #
  check_logical_scalar( trend, "trend")
  
  # -------------------------------------------------------------------------- #
  # 5. Check trend_qua
  # -------------------------------------------------------------------------- #
  check_logical_scalar(trend_qua,"trend_qua")
  

  if (trend_qua && !trend) {
    warning(
      paste0(
        "'trend_qua = TRUE' but 'trend = FALSE'. ",
        "The quadratic trend will be included without a linear trend."
      ),
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------------------------- #
  # 6. Check star
  # -------------------------------------------------------------------------- #
  
  if (!is.null(star)) {
    check_numeric_block(
      x = star,
      argument_name = "star",
      expected_rows = number_of_observations
    )
  }
  
  
  # -------------------------------------------------------------------------- #
  # 7. Check xlags
  # -------------------------------------------------------------------------- #
  check_integer_scalar(xlags, "xlags",minimum = 0L)
  
  # -------------------------------------------------------------------------- #
  # 8. Check det_controls
  # -------------------------------------------------------------------------- #
  if (!is.null(det_controls)) {
    check_numeric_block(
      x = det_controls,
      argument_name = "det_controls",
      expected_rows = number_of_observations
    )
  }
  
  
  # -------------------------------------------------------------------------- #
  # 9. Check ex
  # -------------------------------------------------------------------------- #
  if (!is.null(ex)) {
    check_numeric_block(x = ex,argument_name = "ex",expected_rows = number_of_observations)
  }
  
  # -------------------------------------------------------------------------- #
  # 10. Check ex_lag
  # -------------------------------------------------------------------------- #
  if (!is.null(ex_lag)) {
    check_integer_scalar(ex_lag,"ex_lag",minimum = 0L)
  }
  
  # Interpret ex_lag = NULL as no lags of the
  # exogenous variables.
  if (is.null(ex_lag)) {
    ex_lag <- 0L
  }
  
  # -------------------------------------------------------------------------- #
  # 11. Check feedback_only
  # -------------------------------------------------------------------------- #
  check_logical_scalar(feedback_only,"feedback_only")
  
  # -------------------------------------------------------------------------- #
  # 12. Check consistency across arguments
  # -------------------------------------------------------------------------- #
  
  # xlags has no meaning if star variables were not supplied.
  if (is.null(star) && xlags > 0L) {
    stop(
      paste0("'xlags' is positive, but no star variables were supplied in 'star'."),
      call. = FALSE
    )
  }
  
  # If contemporaneous star variables are excluded,
  # at least one lag must be included.
  if (!is.null(star) && feedback_only && xlags == 0L) {
    stop(
      paste0(
        "'feedback_only = TRUE' requires 'xlags' to be at least 1."
      ),
      call. = FALSE
    )
  }
  
  # feedback_only has no effect when there are no star variables.
  if (is.null(star) && feedback_only) {
    warning(
      paste0("'feedback_only' has no effect because 'star' is NULL."),
      call. = FALSE
    )
  }
  
  # ex_lag has no meaning if exogenous variables were not supplied.
  if (is.null(ex) && ex_lag > 0L) {
    stop(
      paste0("'ex_lag' is positive, but no exogenous variables were supplied in 'ex'."),
      call. = FALSE
    )
  }
  
  # Find the largest lag used anywhere in the model.
  maximum_lag <- max(
    lags,
    if (is.null(star)) 0L else xlags,
    if (is.null(ex)) 0L else ex_lag
  )
  
  # The sample must contain observations after accounting
  # for the largest lag.
  if (number_of_observations <= maximum_lag) {
    stop(
      sprintf(
        paste0(
          "The sample contains %d observations, but the largest lag is %d. ",
          "The number of observations must exceed the largest lag."
        ),
        number_of_observations, maximum_lag
      ),
      call. = FALSE
    )
  }
  
  
  # -------------------------------------------------------------------------- #
  # Normalise optional input blocks
  # -------------------------------------------------------------------------- #
  
  # Convert optional matrices into data frames so that column-name
  # and lag operations behave consistently during matrix construction.
  
  star_data <- if (is.null(star)) {
    NULL
  } else {
    as.data.frame(star,check.names = FALSE)
  }
  
  det_controls_data <- if (is.null(det_controls)) {
    NULL
  } else {
    as.data.frame(det_controls, check.names = FALSE)
  }
  
  ex_data <- if (is.null(ex)) {
    NULL
  } else {
    as.data.frame(ex,check.names = FALSE)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Construction of the model matrices
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  
  # -------------------------------------------------------------------------- #
  # 1. Data preparation and model dimensions
  # -------------------------------------------------------------------------- #
  
  # Convert the endogenous data to a standard data frame.
  y <- as.data.frame(data, check.names = FALSE)
  
  # Number of variables in each original data block.
  number_of_endogenous_variables <- ncol(y)
  
  number_of_star_variables <- if (is.null(star_data)) {
    0L
  } else {
    ncol(star_data)
  }
  
  number_of_det_controls <- if (is.null(det_controls_data)) {
    0L
  } else {
    ncol(det_controls_data)
  }
  
  number_of_exogenous_variables <- if (is.null(ex_data)) {
    0L
  } else {
    ncol(ex_data)
  }
  
  # Number of lagged variables in each block.
  number_of_lagged_endogenous_variables <- number_of_endogenous_variables * lags
  
  number_of_lagged_star_variables <- number_of_star_variables * xlags
  
  number_of_lagged_exogenous_variables <- number_of_exogenous_variables * ex_lag
  
  # Variable names.
  names_of_endog_variables <- names(y)
  
  names_of_star_variables <- if (is.null(star_data)) {
    character(0)
  } else {
    names(star_data)
  }
  
  names_of_det_controls <- if (is.null(det_controls_data)) {
    character(0)
  } else {
    names(det_controls_data)
  }
  
  names_of_exogenous_variables <- if (is.null(ex_data)) {
    character(0)
  } else {
    names(ex_data)
  }
  
  
  # -------------------------------------------------------------------------- #
  # 2. Construct the lagged endogenous-variable block
  # -------------------------------------------------------------------------- #
  
  endogenous_lag_blocks <- lapply(
    seq_len(lags),
    function(lag_order) {
    
      lag_block <- dplyr::lag(y,lag_order)
      names(lag_block) <- paste0(names_of_endog_variables,"_",lag_order)
      
      lag_block
    }
  )
  
  endogenous_lag_block <- do.call(cbind, unname(endogenous_lag_blocks))
  
  number_of_endogenous_lag_regressors <- ncol(endogenous_lag_block)
  
  if (number_of_endogenous_lag_regressors != number_of_lagged_endogenous_variables) {
    stop(
      paste0("Internal error: expected ",number_of_lagged_endogenous_variables,
        " lagged endogenous regressors but constructed ",number_of_endogenous_lag_regressors,"."
      ),
      call. = FALSE
    )
  }
  
  # -------------------------------------------------------------------------- #
  # 3. Construct the deterministic-term block
  # -------------------------------------------------------------------------- #
  time_index <- seq_len(number_of_observations)
  deterministic_blocks <- list()
  # Constant.
  if (const) {
    deterministic_blocks[["constant"]] <- rep(1, number_of_observations)
  }
  # Linear trend.
  if (trend) {
    deterministic_blocks[["trend"]] <- time_index
  }
  # Quadratic trend.
  if (trend_qua) {
    deterministic_blocks[["trend_squared"]] <- time_index^2
  }
  # Additional deterministic controls.
  if (!is.null(det_controls_data)) { 
    deterministic_blocks <- c(deterministic_blocks,as.list(det_controls_data) )
  }
  # Combine all deterministic terms.
  deterministic_block <- if (length(deterministic_blocks) == 0L ) {
    NULL
  } else {
    as.data.frame(deterministic_blocks, check.names = FALSE)
  }
  number_of_deterministic <- if (is.null(deterministic_block) ) {
    0L
  } else {
    ncol(deterministic_block)
  }
  
  # -------------------------------------------------------------------------- #
  # 4. Construct the star-variable block
  # -------------------------------------------------------------------------- #
  star_block <- NULL
  
  if (!is.null(star_data)) {
    
    star_blocks <- list()
    
    # Include contemporaneous star variables unless
    # feedback_only is TRUE.
    if (!feedback_only) {
      
      contemporaneous_star <- star_data
      
      names(contemporaneous_star) <- paste0(names_of_star_variables,"_0")
      
      star_blocks[["contemporaneous"]] <- contemporaneous_star
    }
    
    # Construct the lagged star-variable blocks.
    if (xlags > 0L) {
      
      lagged_star_blocks <- lapply(seq_len(xlags),
                                   
        function(lag_order) {
          
          lag_block <- dplyr::lag(star_data,lag_order)
          
          names(lag_block) <- paste0(names_of_star_variables,"_",lag_order)
          
          lag_block
        }
      )
      
      star_blocks <- c(star_blocks,lagged_star_blocks)
    }
    
    # Combine contemporaneous and lagged star variables.
    star_block <- do.call(cbind, unname(star_blocks))
  }
  
  number_of_contemporaneous_star_variables <-
    if (!is.null(star_data) && !feedback_only) {
      number_of_star_variables
    } else {
      0L
    }
  
  number_of_star_regressors <- if (
    is.null(star_block)
  ) {
    0L
  } else {
    ncol(star_block)
  }
  
  expected_number_of_star_regressors <- number_of_contemporaneous_star_variables + number_of_lagged_star_variables
  
  if (number_of_star_regressors !=  expected_number_of_star_regressors) {
    stop(
      paste0(
        "Internal error: expected ",
        expected_number_of_star_regressors,
        " star regressors but constructed ",
        number_of_star_regressors,
        "."
      ),
      call. = FALSE
    )
  }
  
  
  # -------------------------------------------------------------------------- #
  # 5. Construct the exogenous-variable block
  # -------------------------------------------------------------------------- #
  
  ex_block <- NULL
  
  if (!is.null(ex_data)) {
    
    ex_blocks <- list()
    
    # Contemporaneous exogenous variables.
    contemporaneous_ex <- ex_data
    
    names(contemporaneous_ex) <- paste0(names_of_exogenous_variables, "_0")
    
    ex_blocks[["contemporaneous"]] <- contemporaneous_ex
    
    # Construct the lagged exogenous-variable blocks.
    if (ex_lag > 0L) {
      
      lagged_ex_blocks <- lapply(
        seq_len(ex_lag),
        function(lag_order) {
          
          lag_block <- dplyr::lag(
            ex_data,
            lag_order
          )
          
          names(lag_block) <- paste0(
            names_of_exogenous_variables,
            "_",
            lag_order
          )
          
          lag_block
        }
      )
      
      ex_blocks <- c( ex_blocks, lagged_ex_blocks)
    }
    
    # Combine contemporaneous and lagged exogenous variables.
    ex_block <- do.call(cbind,unname(ex_blocks) )
  }
  
  number_of_exogenous_regressors <- if (
    is.null(ex_block)
  ) {
    0L
  } else {
    ncol(ex_block)
  }
  
  expected_number_of_exogenous_regressors <- number_of_exogenous_variables + number_of_lagged_exogenous_variables
  
  if (
    number_of_exogenous_regressors !=
    expected_number_of_exogenous_regressors
  ) {
    stop(
      paste0(
        "Internal error: expected ",
        expected_number_of_exogenous_regressors,
        " exogenous regressors but constructed ",
        number_of_exogenous_regressors,
        "."
      ),
      call. = FALSE
    )
  }
  
  
  # -------------------------------------------------------------------------- #
  # 6. Combine the right-hand-side blocks
  # -------------------------------------------------------------------------- #
  
  rhs_blocks <- list(
    deterministic = deterministic_block,
    exogenous = ex_block,
    star = star_block,
    endogenous_lags = endogenous_lag_block
  )
  
  # Remove blocks that are NULL.
  rhs_blocks <- rhs_blocks[
    !vapply(
      rhs_blocks,
      is.null,
      logical(1)
    )
  ]
  
  # Combine the remaining blocks without adding
  # the list names to the column names.
  x_full <- do.call(cbind, unname(rhs_blocks))
  
  x_full <- as.data.frame(x_full, check.names = FALSE )
  
  # Ensure that every regressor has a unique name.
  if (anyDuplicated(names(x_full))) {
    
    duplicated_names <- unique(
      names(x_full)[duplicated(names(x_full))]
    )
    
    stop(
      paste0(
        "Duplicated regressor names: ",
        paste(
          duplicated_names,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }
  
  # Number of regressors outside the endogenous-lag block.
  number_of_non_endogenous_regressors <- number_of_deterministic + number_of_exogenous_regressors + number_of_star_regressors
  
  # Total number of regressors per equation.
  number_of_regressors <- ncol(x_full)
  
  expected_number_of_regressors <- number_of_non_endogenous_regressors + number_of_endogenous_lag_regressors
  
  if (
    number_of_regressors !=
    expected_number_of_regressors
  ) {
    stop(
      paste0(
        "Internal error: expected ",
        expected_number_of_regressors,
        " regressors but constructed ",
        number_of_regressors,
        "."
      ),
      call. = FALSE
    )
  }
  
  
  # -------------------------------------------------------------------------- #
  # 7. Define the estimation sample
  # -------------------------------------------------------------------------- #
  
  # Remove the initial observations that are unavailable because of lags.
  estimation_rows <- seq.int(from = maximum_lag + 1L, to = number_of_observations )
  
  # Construct the left-hand-side matrix.
  y_lhs <- as.matrix( y[estimation_rows,  ,drop = FALSE] )
  
  # Construct the right-hand-side matrix.
  x_rhs <- as.matrix( x_full[ estimation_rows,, drop = FALSE]  )
  

  # -------------------------------------------------------------------------- #
  # 8. Check the estimation matrices
  # -------------------------------------------------------------------------- #
  
  # The dependent-variable and regressor matrices must have
  # the same number of observations.
  if (nrow(y_lhs) != nrow(x_rhs)) {
    stop(
      paste0(
        "Internal error: 'y_lhs' and 'x_rhs' ",
        "have different row counts."
      ),
      call. = FALSE
    )
  }
  
  # The number of dependent variables must equal the number
  # of endogenous variables recorded earlier.
  if (
    ncol(y_lhs) !=
    number_of_endogenous_variables
  ) {
    stop(
      paste0(
        "Internal error: expected ",
        number_of_endogenous_variables,
        " endogenous variables but constructed ",
        ncol(y_lhs),
        "."
      ),
      call. = FALSE
    )
  }
  
  # The number of columns in x_rhs must equal the number of
  # regressors recorded earlier.
  if (
    ncol(x_rhs) !=
    number_of_regressors
  ) {
    stop(
      paste0(
        "Internal error: expected ",
        number_of_regressors,
        " regressors but constructed ",
        ncol(x_rhs),
        "."
      ),
      call. = FALSE
    )
  }
  
  # No missing values should remain after removing the
  # initial lagged observations.
  if (anyNA(y_lhs) || anyNA(x_rhs)) {
    stop(
      paste0(
        "Internal error: the estimation matrices contain ",
        "missing values. Check the lag construction and ",
        "'maximum_lag'."
      ),
      call. = FALSE
    )
  }
  
  # All values in the estimation matrices must be finite.
  if (
    !all(is.finite(y_lhs)) ||
    !all(is.finite(x_rhs))
  ) {
    stop(
      paste0(
        "The estimation matrices contain ",
        "non-finite values."
      ),
      call. = FALSE
    )
  }
  
  
  # -------------------------------------------------------------------------- #
  # 9. Sample size and degrees of freedom
  # -------------------------------------------------------------------------- #
  
  number_of_effective_observations <-
    nrow(y_lhs)
  
  dof <-
    number_of_effective_observations -
    number_of_regressors
  
  if (dof <= 0L) {
    warning(
      paste0(
        "The model has ",
        number_of_effective_observations,
        " effective observations and ",
        number_of_regressors,
        " regressors per equation, giving ",
        dof,
        " residual degrees of freedom."
      ),
      call. = FALSE
    )
  }
  
  
  # -------------------------------------------------------------------------- #
  # 10. Align the time identifier
  # -------------------------------------------------------------------------- #
  
  if (!is.null(TimeID)) {
    TimeID <- TimeID[estimation_rows, , drop = FALSE]
    
    if (nrow(TimeID) != number_of_effective_observations) {
      stop(
        paste0("Internal error: 'TimeID' is not aligned with the estimation matrices."),
        call. = FALSE
      )
    }
  }
  
  
  # -------------------------------------------------------------------------- #
  # 11. Return the results
  # -------------------------------------------------------------------------- #
  
  vardata <- list(
    
    # Estimation matrices
    y_lhs = y_lhs,
    x_rhs = x_rhs,
    # Numbers of original variables
    number_of_endogenous_variables =  number_of_endogenous_variables,
    number_of_star_variables = number_of_star_variables,
    number_of_exogenous_variables = number_of_exogenous_variables,
    number_of_det_controls = number_of_det_controls,
    # Numbers of lagged variables
    number_of_lagged_endogenous_variables = number_of_lagged_endogenous_variables,
    number_of_lagged_star_variables = number_of_lagged_star_variables,
    number_of_lagged_exogenous_variables = number_of_lagged_exogenous_variables,
    number_of_contemporaneous_star_variables = number_of_contemporaneous_star_variables,
    # Numbers of regressors
    number_of_deterministic = number_of_deterministic,
    number_of_exogenous_regressors = number_of_exogenous_regressors,
    number_of_star_regressors = number_of_star_regressors,
    number_of_endogenous_lag_regressors = number_of_endogenous_lag_regressors,
    number_of_non_endogenous_regressors = number_of_non_endogenous_regressors,
    number_of_regressors = number_of_regressors,
    # Sizes of the right-hand-side blocks
    rhs_block_sizes = c( 
      deterministic = number_of_deterministic,
      exogenous = number_of_exogenous_regressors,
      star = number_of_star_regressors,
      endogenous_lags = number_of_endogenous_lag_regressors
  ),
    
    # Sample dimensions
    number_of_observations = number_of_observations,
    number_of_effective_observations =number_of_effective_observations,
    estimation_rows = estimation_rows,
    maximum_lag = maximum_lag,
    # Lag specification
    number_of_lags = lags,
    number_of_xlags = xlags,
    number_of_ex_lags = ex_lag,
    feedback_only = feedback_only,
    # Degrees of freedom and dates
    dof = dof,
    TimeID =TimeID,
    # Variable names
    names_of_endog_variables = names_of_endog_variables,
    names_of_star_variables = names_of_star_variables,
    names_of_det_controls = names_of_det_controls,
    names_of_exogenous_variables = names_of_exogenous_variables,
    names_of_rhs_variables = colnames(x_rhs)
  )
  
  return(vardata)
}



