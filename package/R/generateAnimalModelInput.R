#' Generate data input for an animal model.
#'
#' @param formula formula specifying fixed effects and response variables
#' @param data data.frame with data. 
#' @param A relationship matrix from pedigree
#' @param IDcol Which column has the individual IDs which match the A matrix rownames.
#' @param out_type Default is none to output object to R. Choose 'HDF5' or 'matlab' for mat
#' @param out_folder folder to create the Matlab or Julia input files
#' @importFrom R.matlab writeMat
#' @importFrom rhdf5 h5createFile h5createGroup h5write H5close
#' @importFrom stats lm as.formula model.matrix
#' @return
#' list with animal model input objects
#' @export
#' @examples
#' library(ratones)
#' data(ratonesdf)
#' data(ratones_ped)
#' formula = paste0("cbind(",
#'           paste(names(dplyr::select(ratonesdf, IS_PM:BA_OPI)),
#'                 collapse = ", "), ") ~ SEX + AGE")
#' animal_model_data = generateAnimalModelInput(formula, ratonesdf,
#'                                              ratones_ped$A, out_folder = "./")
genAnimalModelInput <- function(formula, data, A, IDcol = "ID", 
                                out_type = c("none", "HDF5", "matlab"),
                                out_folder = "BSFG_run"){
  out_type = match.arg(out_type)
  model = lm(as.formula(formula), data = data)
  Y = model$model[[1]]
  X = model.matrix(model)
  pos = sapply(data[[IDcol]], function(x) which(x == rownames(A)))
  A = as.matrix(A[pos, pos])
  Z = diag(nrow(A))
  if(out_type != "none"){
    if(!dir.exists(out_folder))
      dir.create(file.path(out_folder), showWarnings = FALSE)
    if(out_type == "matlab")
      writeMat(paste0(out_folder, "/setup.mat"), A = A, X = t(X), Y = Y, Z_1 = Z)
    else if(out_type=="HDF5"){
      out_file = file.path(out_folder, "/setup.h5")
      if(file.exists(out_file))
        file.remove(out_file)
      saveHDF5animalModel(Y, X, A, Z, out_file)
    }
  }
  return(list(K = ncol(Y),
              J = ncol(X),
              N = nrow(Y),
              A = A,
              X = X,
              Y = Y,
              Z = Z))
}

#' Generate input for BSFG model in Julia
#'
#' @param Y matrix of observed phenotypes
#' @param X fixed effects model matrix 
#' @param A relationship matrix from pedigree. Use nadiv package.
#' @param Z Which column has the individual IDs which match the A matrix rownames.
#' @param out_file output filename.
#' @importFrom rhdf5 h5createFile h5createGroup h5write H5close
#' @return
#' list with animal model input objects
#' @export
#' @examples
#' library(ratones)
#' saveHDF5animalModel(Y, X, A, Z, "setup.h5")
saveHDF5animalModel = function(Y, X, A, Z, out_file){
  h5createFile(out_file)
  h5createGroup(out_file, "Input")
  h5write(Y, out_file, name="Input/Y")
  h5write(t(X), out_file, name="Input/X")
  h5write(as.matrix(A), out_file, name="Input/A")
  h5write(Z, out_file, name="Input/Z_1")
  H5close()
  return(list(Y, X, A, Z))
}