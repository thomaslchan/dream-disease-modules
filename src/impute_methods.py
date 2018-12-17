import numpy as np


########################### GLOBAL IMPUTATION ##########################
#----------------------------------------------------------------------#


def impute(matrices, imputed_values):
    """
    Replaces zeros by imputing a value per matrix, and refills the diagonals
    with zeros.

    Args:
    -----------------------------------------------------------------
    - matrices: List of matrices as Numpy arrays
    - imputed_values: List of values to impute (same length as matrices)
    """
    for matrix, val in zip(matrices, imputed_values):
        matrix[matrix == 0] = val
        np.fill_diagonal(matrix, 0)
    return matrices


def zero_impute(matrices):
    """
    Simply returns the original aggregated matrices, which are imputed
    with zeroes by default.

    Args:
    -----------------------------------------------------------------
    - matrices: List of distance matrices as Numpy arrays (must be same size)
    """
    return matrices


def mean_impute(matrices):
    """
    Imputes the global mean of each distance matrix and returns the matrices. 

    Args:
    -----------------------------------------------------------------
    - matrices: List of distance matrices as Numpy arrays (must be same size)
    """
    vals = [np.mean(matrix) for matrix in matrices]
    return impute(matrices, vals)


def min_impute(matrices):
    """
    Imputes the global nonzero minimum of each distance matrix 
    and returns the matrices.

    Args:
    -----------------------------------------------------------------
    - matrices: List of distance matrices as Numpy arrays (must be same size)
    """
    vals = [np.min(matrix[np.nonzero(matrix)]) for matrix in matrices]
    return impute(matrices, vals)


def max_impute(matrices):
    """
    Imputes the global max of each distance matrix and returns the matrices.

    Args:
    -----------------------------------------------------------------
    - matrices: List of distance matrices as Numpy arrays (must be same size)
    """
    vals = [np.max(matrix) for matrix in matrices]
    return impute(matrices, vals)


def median_impute(matrices):
    """
    Imputes the global median of each distance matrix and returns the matrices.

    Args:
    -----------------------------------------------------------------
    - matrices: List of distance matrices as Numpy arrays (must be same size)
    """
    vals = [np.median(matrix) for matrix in matrices]
    return impute(matrices, vals)


#============================================================================


######################## LOCAL IMPUTATION ###################################
#----------------------------------------------------------------------------


def preprocess(matrices):
    """
    Flattens list of matrices into 2D matrix with dimensions
    (number of matrices  X  number of elements per matrix). Returns
    the flat matrix and the number of nodes.
    
    Args:
    -----------------------------------------------------------------
    - matrices: List of distance matrices as Numpy arrays
    """
    rows, num_matrices = len(matrices[0]), len(matrices)
    flat = np.vstack(np.dstack(tuple(matrices)))
    flat[flat == 0] = np.nan
    return flat, rows


def mean_local_impute(matrices, global_impute=mean_impute):
    """
    Imputes the elementwise mean across matrices. If there are none,
    we default to using a global impute method. Returns list of matrices.

    Args:
    -----------------------------------------------------------------
    - matrices: List of matrices
    - global_impute: Global impute function (default=mean_impute)
    """
    flat, rows = preprocess(matrices)
    mean = np.nanmean(flat, axis=1).reshape(rows, rows)
    vals = global_impute([np.nan_to_num(mean)])
    return local_impute(matrices, vals)


def median_local_impute(matrices, global_impute=mean_impute):
    """
    Imputes the elementwise median across matrices. If there are none,
    we default to using a global impute method. Returns list of matrices.

    Args:
    -----------------------------------------------------------------
    - matrices: List of matrices
    - global_impute: Global impute function (default=mean_impute)
    """
    flat, rows = preprocess(matrices)
    mean = np.nanmedian(flat, axis=1).reshape(rows, rows)
    vals = global_impute([np.nan_to_num(mean)])
    return local_impute(matrices, vals)


def min_local_impute(matrices, global_impute=mean_impute):
    """
    Imputes the elementwise min across matrices. If there are none,
    we default to using a global impute method. Returns list of matrices.

    Args:
    -----------------------------------------------------------------
    - matrices: List of matrices
    - global_impute: Global impute function (default=mean_impute)
    """
    flat, rows = preprocess(matrices)
    mean = np.nanmin(flat, axis=1).reshape(rows, rows)
    vals = global_impute([np.nan_to_num(mean)])
    return local_impute(matrices, vals)


def max_local_impute(matrices, global_impute=mean_impute):
    """
    Imputes the elementwise max across matrices. If there are none,
    we default to using a global impute method. Returns list of matrices.

    Args:
    -----------------------------------------------------------------
    - matrices: List of matrices
    - global_impute: Global impute function (default=mean_impute)
    """
    flat, rows = preprocess(matrices)
    mean = np.nanmax(flat, axis=1).reshape(rows, rows)
    vals = global_impute([np.nan_to_num(mean)])
    return local_impute(matrices, vals)


def local_impute(matrices, vals):
    """
    Imputes some matrix of values for each corresponding missing element
    in the matrices. Returns the list of matrices.

    Args:
    -----------------------------------------------------------------
    - matrices: List of matrices
    - global_impute: Global impute function (default=mean_impute)
    """
    [np.copyto(matrix, vals, where=(matrix==0)) for matrix in matrices]
    [np.fill_diagonal(matrix, 0) for matrix in matrices]
    return matrices

#====================================================================