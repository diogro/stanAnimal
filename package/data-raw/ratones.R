# Builds the bundled mouse datasets (data/ratones_morpho.rda and
# data/ratones_pedigree.rda) from the companion data package diogro/ratones.
#
# The original records come from Penna et al. (2017), archived at Dryad
# (doi:10.5061/dryad.5gr8r) and packaged in https://github.com/diogro/ratones,
# which is MIT licensed. Set RATONES_SOURCE to the data/ directory of a
# checkout of that repository, or leave it unset to use the default location.
#
# Run from the package directory: Rscript data-raw/ratones.R

source_dir <- Sys.getenv("RATONES_SOURCE",
                         unset = file.path("..", "..", "ratones", "data"))
load(file.path(source_dir, "ratonesdf.RData"))
load(file.path(source_dir, "ratones_ped.RData"))

ratones_morpho <- ratonesdf
ratones_morpho$ID <- as.integer(ratones_morpho$ID)
ratones_morpho$GER <- as.integer(ratones_morpho$GER)
ratones_morpho$MADRE <- as.integer(ratones_morpho$MADRE)
ratones_morpho$PADRE <- as.integer(ratones_morpho$PADRE)
# Row names carry the individual IDs, so animal_design() can match a pedigree
# or a relationship matrix to these rows without any extra preparation.
rownames(ratones_morpho) <- as.character(ratones_morpho$ID)

ped <- ratones_ped$ped
ratones_pedigree <- data.frame(id = as.integer(ped$ID),
                               sire = as.integer(ped$Sire),
                               dam = as.integer(ped$Dam))

stopifnot(
  identical(dim(ratones_morpho), c(329L, 46L)),
  identical(dim(ratones_pedigree), c(525L, 3L)),
  !anyDuplicated(ratones_morpho$ID),
  !anyDuplicated(ratones_pedigree$id),
  all(ratones_morpho$ID %in% ratones_pedigree$id)
)

dir.create("data", showWarnings = FALSE)
save(ratones_morpho, file = "data/ratones_morpho.rda", compress = "xz", version = 2)
save(ratones_pedigree, file = "data/ratones_pedigree.rda", compress = "xz", version = 2)
