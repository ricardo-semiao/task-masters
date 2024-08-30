# Tasks From EESP-FGV's Masters in Economics

Welcome! This is the repository for the tasks I've done during my Masters in Economics at EESP-FGV. Each folder belongs to a different class, but have a common structure:

- Each problem set/task has its own folder, that all contain similar files:
    - A _psi\_instructions.pdf_ file, with the questions of the P.S.
    - A _psi.R_ and _psi.tex_ files, to generate the results and the final report, _psi.pdf_.
    - A _data_, _figures_, and _tables_ folders, with assets.
    - A _latex\_build_ folder, with auxiliary files for the .tex file compilation.

Below, I talk about each of the classes' tasks.

## [Econometrics II](econometrics2)

An time series econometrics class. The problems were varied, with math derivations, empirical analysis, and simulations. The themes were:

- Stationarity analysis.
- ARIMA, VAR, VECM models.
- GMM.
- Dynamic Panel Models.


## [Forecasting](forecasting)

Also an time series econometrics task, but focused on forecasting. Here i've used other techniques such as:

- Recursive OLS.
- Model/forecasts evaluation and combination.
- Model diagnostics.
- High dimensionality models, Bayesian VARs.


## [Quantitative Methods in Macroeconomics](quant-macro)

A class on numeric methods to solve macroeconomic models. The tasks were about:

- Basics of numeric methods.
- Methods such as Bisection, Newton-Rapson, Secant, etc.
- Solving dynamic programming problems and Value Function Iteration.
- Macro specific problems, such as solving a RBC model.


## [Econometrics I](econometrics1)

A introductory class on econometrics. The task is described at [instructions.pdf](econometrics1/instructions.pdf), but basically entails the replication of works in the literature of “quantity-quality trade-off”, a la Becker and Lewis (1973) and Becker and Tomes (1976). The finished report is [text.pdf](econometrics1/text.pdf).

The data comes from IPUMS International, the Puerto Rico 2010's census, and is at [data/](econometrics1/data/). The main script is [main.R](econometrics1/main.R), and the report is built with [text.tex](econometrics1/text.tex).



## [Microeconometrics I](microeconometrics1)

A class focused on experiments and causal inference. Again, very varied, mixing theoretical and empirical problems. The tasks were about:

- Identification.
- Experiments attrition, balance, and diagnostics.
