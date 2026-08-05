Currently under development.

Installation of dependencies currently needs to be done by hand. To use the emulation options requires at least 1 of the following:

```{r}
install.packages(RobustGaSP)
install.packages(hetGP)
devtools::install_github("mingdeyu/dgpsi-R")
```

You may also need to explicitly load:

```{r}
library(ggplot2)
```

The examples in `tutorials/` should work, but some of the functions may change slightly at a later date.

