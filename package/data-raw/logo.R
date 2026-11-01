# Hex-sticker logo for stanimal.
#
# Cuts the artwork in the repository root (stanimal.png) to its hexagonal
# outline and writes man/figures/logo.png with a transparent background -- the
# file pkgdown picks up as the package logo and that the README shows.
#
# Run from the package directory:
#   Rscript data-raw/logo.R
#
# Needs the png package.

repo <- normalizePath(file.path(normalizePath("."), ".."))
input <- file.path(repo, "stanimal.png")
output <- file.path("man", "figures", "logo.png")

if (!file.exists(input))
  stop("missing ", input, ": convert stanimal.jpeg to PNG first")
if (!requireNamespace("png", quietly = TRUE))
  stop("the png package is needed to write the logo")

art <- png::readPNG(input)
height <- dim(art)[1]
width <- dim(art)[2]

# The drawing is a hexagon on a white canvas, so its extent is the extent of
# the non-white (outline) pixels.
ink <- art[, , 1] < 0.97 | art[, , 2] < 0.97 | art[, , 3] < 0.97
rows <- which(rowSums(ink) > 0)
cols <- which(colSums(ink) > 0)
ys <- range(rows)
xs <- range(cols)

# Regular hexagon, one vertex at the top and one at the bottom, inscribed in
# that extent.
centre_x <- mean(xs)
centre_y <- mean(ys)
half_width <- diff(xs) / 2
half_height <- diff(ys) / 2
hexagon <- cbind(
  x = c(centre_x, centre_x + half_width, centre_x + half_width,
        centre_x, centre_x - half_width, centre_x - half_width),
  y = c(centre_y - half_height, centre_y - half_height / 2,
        centre_y + half_height / 2, centre_y + half_height,
        centre_y + half_height / 2, centre_y - half_height / 2)
)

# Signed distance from each pixel centre to the hexagon boundary, positive
# outside. Everything outside is cut; ramping alpha across one pixel keeps the
# cut edge smooth. The pixels are laid out rows first, as in the image array.
points <- expand.grid(y = seq_len(height), x = seq_len(width))
distance <- matrix(-Inf, height, width)
for (i in seq_len(nrow(hexagon))) {
  from <- hexagon[i, ]
  to <- hexagon[if (i == nrow(hexagon)) 1L else i + 1L, ]
  edge <- to - from
  normal <- c(edge[2], -edge[1])
  normal <- normal / sqrt(sum(normal^2))
  # Point the normal away from the centre of the hexagon.
  if (sum((c(centre_x, centre_y) - from) * normal) > 0) normal <- -normal
  d <- (points$x - from[1]) * normal[1] + (points$y - from[2]) * normal[2]
  distance <- pmax(distance, matrix(d, height, width))
}

# The artwork fades from the dark outline to white over a few pixels, and its
# hexagon has rounded corners that leave white wedges just inside a sharp
# hexagon. Cut well inside the outer edge so that the last visible pixels are
# the outline itself, instead of a pale rim that shows up on dark backgrounds.
# The rim scales with the artwork, hence the fraction of its width.
trim <- 0.014 * width

logo <- array(1, c(height, width, 4))
logo[, , 1:3] <- art[, , 1:3]
logo[, , 4] <- pmin(pmax(0.5 - (distance + trim), 0), 1)

logo <- logo[ys[1]:ys[2], xs[1]:xs[2], , drop = FALSE]

# The artwork is drawn far larger than the README or the pkgdown navbar ever
# show it, so halve it once. Averaging has to happen on premultiplied colours:
# the transparent pixels still carry white, and averaging that in would put a
# pale rim back around the sticker.
premultiply <- function(x) {
  for (k in 1:3) x[, , k] <- x[, , k] * x[, , 4]
  x
}
unpremultiply <- function(x) {
  alpha <- x[, , 4]
  solid <- alpha > 0
  for (k in 1:3) x[, , k][solid] <- pmin(1, x[, , k][solid] / alpha[solid])
  x
}
halve <- function(x) {
  rows <- 2 * (dim(x)[1] %/% 2)
  cols <- 2 * (dim(x)[2] %/% 2)
  first <- seq(1, rows, 2)
  second <- seq(2, rows, 2)
  left <- seq(1, cols, 2)
  right <- seq(2, cols, 2)
  (x[first, left, , drop = FALSE] + x[second, left, , drop = FALSE] +
   x[first, right, , drop = FALSE] + x[second, right, , drop = FALSE]) / 4
}

logo <- unpremultiply(halve(premultiply(logo)))

dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
png::writePNG(logo, output)
cat("wrote", output, "-", dim(logo)[2], "x", dim(logo)[1], "pixels\n")
