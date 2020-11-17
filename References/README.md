# Relevant methods:

We will compare the following methods using the US spatio-temporal data:

- PROPOSAL: Our model, with a set of covariates and with spatio-temporal covariance structure defined by exponentially decaying correlation.
- BASE: Linear model with only intercept term.
- LM: Linear model with covariates same as our proposed model.
- ARIMA: Define appropriate ARIMA models (chosen by `auto.arima` function) for every location separately.
- SEIR: Generalized SEIR model, refer to [the paper by Peng et al](https://github.com/soudeepd/COVID_Analysis/blob/main/References/Epidemic%20analysis%20of%20COVID-19%20in%20China%20bydynamical%20modeling.pdf).
- SLX: Spatial Durbin linear model with lagged regressors. Refer to [this link for implementation](https://r-spatial.github.io/spatialreg/reference/SLX.html#examples) and [the paper by Guliyev](https://github.com/soudeepd/COVID_Analysis/blob/main/References/Determining%20the%20spatial%20effects%20of%20COVID-19%20using%20the%20spatial%20panel%20data%20model.pdf) for related results and discussion.
