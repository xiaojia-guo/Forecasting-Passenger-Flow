## Code for Forecasting Transfer Passenger Flows at Heathrow

This code is used in our paper "Forecasting Airport Passenger Flows Using Real-Time Data and Machine Learning". In collaboration with Heathrow airport, we develop a predictive system that generates quantile forecasts of transferring passengers' connection times, and the number of passengers arriving at the immigrating and security areas.

The predictive model developed is based on a regression tree combined with copula-based simulations. We generalize the tree method to predict complete distributions, moving beyond point forecasts.

There are xxx files in this repository: 
\\
Fit regression tree.ipynb  ---  Cross-validation for choosing the tuning parameters, fit the tree to the entire training set,     
                                visualize the reduced version of the tree, and find the segment that each passenger in the 
                                test set belongs to.
Pax_flow_simulation.R      ---  Sample connection times from the gamma distributions, and generate quantile forecasts for the 
                                passengers flows, or the number of arrivals in every 15-minute interval. This code als 
                                provides the accuracy test of the connection times and the passenger flows in the test set.
