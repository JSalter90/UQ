# UQ

Code for different emulation/calibration work

* `code/` contains general emulation/calibration code used across multiple projects

  * `code/FastHM.R` contains fast history matching code (for high dimensional output) from Salter \& and Williamson 2022 https://www.dl.begellhouse.com/journals/52034eb04b657aea,7ba373b25cc40b46,0c41fc6d5a39a9c2.html
  
  * `code/rotation_functions` contains main basis emulation code from Salter et al 2019 https://www.tandfonline.com/doi/shareview/10.1080/01621459.2018.1514306 Examples relating to this code are also in this folder, as well as wrappers for combining this code with RobustGasp and mogp_emulator (reproduced from old repository https://github.com/JSalter90/RotateBasis)

* `papers/` contains specific applications (usually related to a paper)

See also mogp_emulator, BayesExeter, ExeterUQ.








