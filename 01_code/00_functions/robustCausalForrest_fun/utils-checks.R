#' Check that a value is a positive integer
#'
#' Internal validation helper.
#'
#' @param x Object to check.
#' @param name Name of the argument (used in error message).
#' @noRd
#'

.require_arg <- function(x, name, context = NULL) {
    if (is.null(x)) {
        msg <- if (is.null(context)) {
            sprintf("`%s` must be specified.", name)
        } else {
            sprintf("For %s, `%s` must be specified.", context, name)
        }

        abort(
            msg,
            class = "rcf_missing_argument",
            arg = name
        )
    }
}

.check_positive_integer <- function(x, name) {

    if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
        rlang::abort(
            message = sprintf("`%s` must be a single positive integer.", name),
            class = "rcf_invalid_argument",
            arg = name,
            value = x
        )
    }

    if (x != as.integer(x)) {
        rlang::abort(
            message = sprintf("`%s` must be an integer.", name),
            class = "rcf_invalid_argument",
            arg = name,
            value = x
        )
    }
}

.check_ratio <- function(x, name) {

    if (!is.numeric(x) || length(x) != 1 || !is.finite(x) ||
        x <= 0 || x > 1) {

        rlang::abort(
            message = sprintf("`%s` must be a single numeric value in (0, 1].", name),
            class = "rcf_invalid_argument",
            arg = name,
            value = x
        )
    }
}
