Currently under development.

The examples in `tutorials/` will run, but some of the functions may be improved, and hence functionality/options change slightly, at a later date.

As shown in the examples, each of the .R source files needs calling:

```{r}
source('plot.R')
source('basis.R')
source('emulate.R')
source('calibrate.R')
```

Installation of dependencies currently needs to be done by hand. To use the emulation options requires at least 1 of the following:

```{r}
install.packages("RobustGaSP")
install.packages("hetGP")
devtools::install_github("mingdeyu/dgpsi-R")
```

You will also need to explicitly load:

```{r}
library(ggplot2)
```


