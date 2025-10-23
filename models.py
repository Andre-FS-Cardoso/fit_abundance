import statsmodels.api as sm
import piecewise_regression
import numpy as np

def fit_models(x_array, y_array, ey_array):

    x = np.array(x_array)
    y = np.array(y_array)
    ey = np.array(ey_array)
    
    # --- REMOVE NaNs ---
    mask = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(ey)
    x = x[mask]
    y = y[mask]
    ey = ey[mask]
    
    if len(x) >= 10:
    
        # CASE 1 fit: simple linear regression
        X = sm.add_constant(x)
        model = sm.OLS(y, X)
        results = model.fit()
        a2 = results.params[1]
        ea2 = results.bse[1]
        b0 = results.params[0]
        eb0 = results.bse[0]
        RSS1 = results.ssr
        
        # CASE 2 fit: 1 breakpoint
        fit2 = piecewise_regression.main.Fit(x, y, n_boot=200, n_breakpoints=1, min_distance_to_edge=0.05)
        RSS2 = fit2.get_results()["rss"] if fit2.get_results()["converged"] else 1e6
        
        # CASE 3 fit: 2 breakpoints
        fit3 = piecewise_regression.main.Fit(x, y, n_boot=200, n_breakpoints=2, 
                                             min_distance_between_breakpoints=0.20, min_distance_to_edge=0.05,
                                             start_values=[0.5, 1.5])
        RSS3 = fit3.get_results()["rss"] if fit3.get_results()["converged"] else 1e6

        # Functions for AIC
        def llf_(X, rss):
            nobs = float(X.shape[0])
            llf = -0.5 * nobs * (np.log(2*np.pi) + np.log(rss/nobs) + 1)
            return llf
            
        def aic_(X, rss, param):
            llf = llf_(X, rss)
            return -2*llf + 2*param
        
        def aic_c(X, rss, param):
            llf = llf_(X, rss)
            return -2*llf + 2*param + (2*param*(param+1))/(len(X)-param-1)
            
        def aic_final(X, RSS, param):
            nobs = float(X.shape[0])
            if nobs/param >= 40:
                AIC = aic_(X, RSS, param)
            else:
                AIC = aic_c(X, RSS, param)
            return AIC
            
        AIC1 = aic_final(x, RSS1, 2)
        AIC2 = aic_final(x, RSS2, 4)
        AIC3 = aic_final(x, RSS3, 6)

        # Selection of the best model
        AICs = [AIC1, AIC2, AIC3]
        best_case = np.argmin(AICs) + 1

        return {
            'x': x, 'y': y, 'ey': ey,
            'best_case': best_case,
            'fit1': (a2, ea2, b0, eb0, RSS1),
            'fit2': fit2,
            'fit3': fit3,
            'AICs': AICs
        }

