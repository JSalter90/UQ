## Getting started

Currently not set up as a real `R` package, so a little work is required to set things up.

Installation of dependencies currently needs to be done by hand. To use the emulation options requires at least 1 of the following:

```{r}
install.packages("RobustGaSP")
install.packages("hetGP")
devtools::install_github("mingdeyu/dgpsi-R")
```

The default emulation code uses `RobustGaSP`, so installation of this is recommended. Others are optional.

To use plotting functionality requires loading `ggplot2`:

```{r}
library(ggplot2)
```

There is a set of examples in `tutorials/`. All should run fine, although some of the functionality/options may change slightly at a later date.

To get started, as shown in each of the examples, each of the .R source files needs calling:

```{r}
source('plot.R')
source('basis.R')
source('emulate.R')
source('calibrate.R')
```
