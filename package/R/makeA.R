#' Additive relationship matrix from a pedigree
#'
#' Computes the numerator (additive) relationship matrix \code{A} from a
#' pedigree, for use as the \code{A} argument of \code{\link{lmm_animal}}.
#'
#' @param pedigree Data frame whose first three columns are the individual, its
#'   sire and its dam, with \code{0}, \code{NA} or \code{""} for unknown
#'   parents. Individuals can be in any order, and the pedigree can contain
#'   ancestors that have no records.
#' @return A symmetric positive-definite matrix with the individuals' IDs on both
#'   dimensions, in the order of \code{pedigree}.
#' @details The additive relationship is computed with the tabular method: with
#'   \eqn{s_i} and \eqn{d_i} the parents of individual \eqn{i},
#'   \deqn{A_{ii} = 1 + A_{s_i d_i} / 2, \qquad A_{ij} = (A_{s_i j} + A_{d_i j}) / 2,}
#'   taking a term as zero when the parent is unknown. The first line includes
#'   the inbreeding coefficient of \eqn{i} through the relationship between its
#'   parents, so pedigrees with inbreeding are handled. Individuals are
#'   processed in an order where parents come before their offspring, whatever
#'   order the pedigree is in, and a pedigree where an individual is its own
#'   ancestor is an error.
#'
#'   A parent that is not one of the individuals of the pedigree is added to
#'   \code{A} as a founder, with a warning, so a truncated pedigree or a
#'   mistyped ID does not silently drop the relationship. Subset \code{A} by
#'   name to match the rows of \code{Y}, or let \code{\link{animal_design}} do
#'   it.
#'
#'   \code{\link{animal_design}} calls this function when it is given
#'   \code{pedigree} and no genotypes.
#' @seealso \code{\link{genomic_relationship}}, \code{\link{animal_design}}
#' @examples
#' ped <- data.frame(id = 1:6, sire = c(0, 0, 1, 1, 4, 4),
#'                   dam = c(0, 0, 2, 2, 3, 5))
#' A <- makeA(ped)
#' round(A, 3)
#' # full-sib mating: individual 6 is the offspring of the sibs 4 and 5
#' ped2 <- data.frame(id = 1:6, sire = c(NA, NA, NA, 1, 1, 4),
#'                    dam = c(NA, NA, NA, 2, 2, 5))
#' diag(makeA(ped2))   # 1.25 for the inbred individual
#' @export
makeA <- function(pedigree) {
  ped <- .clean_pedigree(pedigree)
  ids <- as.character(ped$id)
  if (anyDuplicated(ids))
    stop("`pedigree` has duplicated individuals.", call. = FALSE)

  parent_ids <- unique(c(as.character(ped$sire), as.character(ped$dam)))
  parent_ids <- parent_ids[!is.na(parent_ids)]
  added <- setdiff(parent_ids, ids)
  if (length(added))
    warning("parents that are not individuals in the pedigree were added as ",
            "founders: ", paste(utils::head(added, 5L), collapse = ", "),
            if (length(added) > 5L) ", ..." else "", call. = FALSE)

  all_ids <- c(ids, added)
  n <- length(all_ids)
  # founders added for missing parents have no parents of their own
  founder <- rep(0L, length(added))
  sire <- c(.parent_positions(ped$sire, all_ids), founder)
  dam <- c(.parent_positions(ped$dam, all_ids), founder)
  order <- .pedigree_order(sire, dam, n)

  A <- diag(n)
  processed <- integer(n)
  count <- 0L
  for (position in seq_len(n)) {
    i <- order[position]
    if (count) {
      previous <- processed[seq_len(count)]
      from_sire <- if (sire[i] > 0L) A[sire[i], previous] else 0
      from_dam <- if (dam[i] > 0L) A[dam[i], previous] else 0
      values <- 0.5 * (from_sire + from_dam)
      A[i, previous] <- values
      A[previous, i] <- values
    }
    if (sire[i] > 0L && dam[i] > 0L)
      A[i, i] <- 1 + 0.5 * A[sire[i], dam[i]]
    count <- count + 1L
    processed[count] <- i
  }

  dimnames(A) <- list(all_ids, all_ids)
  A
}

# The first three columns of a pedigree, named id, sire and dam, with unknown
# parents coded as NA.
.clean_pedigree <- function(pedigree) {
  ped <- as.data.frame(pedigree, stringsAsFactors = FALSE)
  if (ncol(ped) < 3L)
    stop("`pedigree` must have at least three columns: individual, sire and dam.",
         call. = FALSE)
  ped <- ped[, 1:3, drop = FALSE]
  names(ped) <- c("id", "sire", "dam")
  ped[] <- lapply(ped, function(x) if (is.factor(x)) as.character(x) else x)
  if (anyNA(ped$id))
    stop("`pedigree` has individuals with a missing ID.", call. = FALSE)
  for (parent in c("sire", "dam")) {
    x <- ped[[parent]]
    unknown <- if (is.numeric(x)) !is.na(x) & x == 0
               else !is.na(x) & x %in% c("", "0")
    x[unknown] <- NA
    ped[[parent]] <- x
  }
  ped
}

# Positions of the parents in the ID vector, with 0 for an unknown parent.
.parent_positions <- function(x, ids) {
  position <- match(as.character(x), ids)
  position[is.na(position)] <- 0L
  position
}

# Order in which individuals can be processed so that parents come before their
# offspring.
.pedigree_order <- function(sire, dam, n) {
  n_parents <- (sire > 0L) + (dam > 0L)
  children <- vector("list", n)
  for (i in seq_len(n)) {
    for (parent in c(sire[i], dam[i])) {
      if (parent > 0L) children[[parent]] <- c(children[[parent]], i)
    }
  }

  ready <- which(n_parents == 0L)
  order <- integer(n)
  filled <- 0L
  head <- 1L
  while (head <= length(ready)) {
    i <- ready[head]
    head <- head + 1L
    filled <- filled + 1L
    order[filled] <- i
    for (child in children[[i]]) {
      n_parents[child] <- n_parents[child] - 1L
      if (n_parents[child] == 0L) ready <- c(ready, child)
    }
  }
  if (filled < n)
    stop("`pedigree` has a loop: some individuals appear as their own ",
         "ancestors.", call. = FALSE)
  order
}
