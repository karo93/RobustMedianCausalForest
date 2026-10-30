#' Control function for inference options for BLB and IJ for the \code{robustCausalForest}
#'
#' @param varianceMethod Character vector specifying the variance estimation
#'   method(s) to use. Possible values are `"IJ"` (Infinitesimal Jackknife)
#'   and `"BLB"` (Bag of Little Bootstraps). Defaults to both.
#' @param ij A named list of additional arguments passed to the
#'   Infinitesimal Jackknife (IJ) variance procedure. Defaults to an empty list.
#' @param blb A named list of arguments controlling the
#'   Bag of Little Bootstraps (BLB) procedure. Elements include:
#'   \describe{
#'     \item{ntrees.bag}{Integer. Total number of trees grown for the BLB method.
#'       Default is 5000.}
#'     \item{nbag}{Integer. Number of little bags (subsamples). Default is 10.}
#'     \item{s_ratio}{Numeric in (0,1]. Subsample size ratio used to construct
#'       each little bag. Default is 0.8.}
#'     \item{minsize_blb}{Integer. Minimum leaf size used when fitting trees
#'       within each BLB subsample. Default is 2.}
#'   }
#' @param ... Additional arguments passed to downstream inference procedures.
#'
#' @return An object of class `"rcf_inference_options"`, a named list containing
#'   the selected variance method(s) (`varianceMethod`) and configuration lists
#'   for the IJ (`ij`) and BLB (`blb`) procedures.
#'
#' @export
#'

inference_options <- function(
        varianceMethod = c("IJ", "BLB"),
        ij = list(),
        blb = list(
            ntrees.bag = 4e3,
            nbag = 10,
            s_ratio = 0.8,
            minsize_blb = 5
        ),
        ...) {

    varianceMethod <- rlang::arg_match(
        varianceMethod,
        values = c("IJ", "BLB"),
        multiple = TRUE
    )

    # ---- BLB checks ----
    if ("BLB" %in% varianceMethod) {
        .require_arg(blb$nbag, "nbag", "varianceMethod = \"BLB\"")
        .require_arg(blb$s_ratio, "s_ratio", "varianceMethod = \"BLB\"")

        .check_positive_integer(blb$nbag, "nbag")
        .check_ratio(blb$s_ratio, "s_ratio")

        if (!is.null(blb$ntrees.bag))
            .check_positive_integer(blb$ntrees.bag, "ntrees.bag")

        if (!is.null(blb$minsize_blb))
            .check_positive_integer(blb$minsize_blb, "minsize_blb")
    }

    structure(
        list(
            varianceMethod = varianceMethod,
            ij = ij,
            blb = blb
        ),
        class = "rcf_inference_options"
    )
}
