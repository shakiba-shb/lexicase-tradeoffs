# lexicase-tradeoffs
This repository contains code and supplimental material for "Theoretical Limits on the Success of Lexicase Selection Under Contradictory Objectives" which is submitted at <a href="https://gecco-2023.sigevo.org/HomePage">The Genetic and Evolutionary Computation Conference</a>, 2023. 


  **Contents:**
  - [Abstract](https://github.com/shakiba-shb/lexicase-tradeoffs#abstract)
  - [Authors](https://github.com/shakiba-shb/lexicase-tradeoffs#authors)
  - [Dependencies](https://github.com/shakiba-shb/lexicase-tradeoffs#dependencies)
  - [Compilation](https://github.com/shakiba-shb/lexicase-tradeoffs#compilation)

  
## Abstract
Lexicase selection is a state of the art parent selection technique for problems that can be broken down into multiple selection criteria. Originally, lexicase selection was developed for cases where these selection criteria are not likely to be in conflict with each other, but it is unclear how critical the assumption that selection criteria do not contradict is. Some prior work has found cases where lexicase selection fails to find a Pareto-optimal solution due to the presence of multiple objectives that contradict each other, while in other cases, lexicase selection has performed well despite the presence of such objectives. Here, we develop theory identifying circumstances under which lexicase selection will or will not fail to find a Pareto-optimal solution. To make this problem tractable, we restrict our investigation to a theoretical problem with maximal contradiction in its objectives. Ultimately, we find that lexicase selection can perform well under many circumstances involving contradictory objectives, but that there are limits to the parameter spaces where high performance is possible. Additionally, we find that epsilon-lexicase selection is much more strongly impacted by contradictory objectives. Our results inform parameter value decisions under lexicase selection and decisions about which problems to use lexicase selection for.

## Authors
- [Shakiba Shahbandegan](https://github.com/shakiba-shb)
- [Emily Dolson](http://emilyldolson.com/)

## Compilation
The results of this project are reported in 2 parts: math results and experimental results.

### Math results:
To simulate our math results, you can use the "lexicase_tradeoffs.py" file. You must install the python library pybind11 to make sure you can use functions adopted from C++ code. For more information about pybind11 installation and use click <a href="https://pybind11.readthedocs.io/en/stable/installing.html"> here.</a> 

```{bash, shell_installation}

# Clone the repo for this project
git clone --recursive https://github.com/shakiba-shb/lexicase-tradeoffs.git

# Compile EC ecology toolbox
cd ec_ecology_toolbox
make

# Set your selected parameters and run the main code file:
python lexicase_tradeoffs.py

```
### Experimental results:
We used the Modular Agent Based Evolver (MABE) software framework for our experimental results. Specifically we used MABE version 2.0 which provides the ability to construct diverse computational evolution scenarios. For more information on how to use MABE 2.0 click <a href="https://mabe2.readthedocs.io/en/latest"> here.</a>

## Results:
- Whether lexicase selection succeeds or fails in solving a multi-objective optimization problem depends on the ratio of genotype dimension and the population size. If these values do not satisfy the inequality in Equation (5), leixcase selection fails to find any optima. 
- Our $\epsilon$-lexicase selection results suggest that higher values of $\epsilon$ reduce lexicase selection's ability to find Pareto-optimal solutions, even when the population size is large relative to dimension. Moreover, this effect is much more dramatic when the genome is composed of floating point values than when it is composed of integers.
- Increasing the value of $\epsilon$ further impedes $\epsilon$-lexicase selection's ability to find Pareto-optimal solutions.




