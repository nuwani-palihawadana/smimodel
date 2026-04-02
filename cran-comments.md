## Resubmission

This is a resubmission following feedback from CRAN.

* Examples and unit tests using the Gurobi solver now run conditionally.
* gurobi is listed under Suggests in the DESCRIPTION.
* Explained how to access the Gurobi solver in the Description field of DESCRIPTION, since the gurobi R package is not available in a CRAN-like repository.

Additionally,
* The verbose argument in SMI model functions has been changed to a named list to provide more flexible control over verbosity options.