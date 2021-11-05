## Resubmission
This is a resubmission. In this version I have:

* I have inserted the references of the used methods in the DESCRIPTION file in the form authors (year) <doi:...>
* In the example section I have replaced \dontrun with \donttest. Unfortunately, there is no way to obtain a execution time less than 5 seconds. The resampling takes time and the data set is with a N of 47 already small. Even the examples with only 250 bootstrap repetitions took more than 5 seconds. Thank you for your understanding.
* The option (warn=-1) is removed and the code revised so that no warning occures any more.
* I do not longer install any packages in my functions, there are paste to the Imports argument in the DESCRIPTION
* I revised the code and it is now possible for the user to manually select the used number of cores. I did that for the reason, that the examples are run on only two cores.
* Concerning possible misspellings: al and et are part of "et al." for referencing more than two authors.

Thanks you for your patience to check my resubmission. Best wishes Sebastian


## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
