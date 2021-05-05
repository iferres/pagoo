# pagoo 0.3.9

* Added panaroo_2_pagoo function to read the output of the panaroo pangenome reconstruction software.

# pagoo 0.3.8

* Improve backward compatibility. Older pagoo objects created by third party packages which depend on pagoo do not have an attribute required to successfully load them into session. Now the approach is to downgrade the object to a base pagoo class, or to provide the namespace using the pkg argument.

# pagoo 0.3.7

* Fixed missing link in documentation.
* Smaller toy dataset to comply with CRAN policies.

# pagoo 0.3.6

* Bugfix, adds compatibility between previous save/load_pangenomeRDS methods

# pagoo 0.3.5

* Faster roary_2_pagoo()
* Improved save_pangenomeRDS() and load_pangenomeRDS() methods
* Several bugfixes.

# pagoo 0.3.3.9000

* Added a `NEWS.md` file to track changes to the package.
