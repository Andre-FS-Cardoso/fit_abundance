import numpy as np
import pandas as pd
import os

def points(name, criterion, r, OH, OH_err, EWHa, Ha6562_cor, OIII5006_cor, NII6583_cor, calibrator, save_table):
    """
    Applies spectral filters based on different criteria and returns the numpy arrays x (radius), y (OH.O3N2.PP04), and yerr (error in OH).
    """

    # --- Mask of valid values ---
    mask_valid = (
        np.isfinite(OH) &
        np.isfinite(EWHa) &
        np.isfinite(OH_err)
    )

    # Apply initial mask to all arrays
    OH       = OH[mask_valid]
    OH_err   = OH_err[mask_valid]
    EWHa     = EWHa[mask_valid]
    Ha6562_cor = Ha6562_cor[mask_valid]
    NII6583_cor = NII6583_cor[mask_valid]
    OIII5006_cor = OIII5006_cor[mask_valid]
    r = r[mask_valid]

    # --- Calculate useful ratios ---
    X_N2Ha = NII6583_cor - Ha6562_cor
    OIII5006 = OIII5006_cor
    EW = EWHa

    # --- Specific filters by criterion ---
    if criterion is None or criterion.lower() == 'none':
        mask = np.ones_like(OH, dtype=bool)
    elif criterion == 'ST06':
        mask = (X_N2Ha <= -0.30) & (OIII5006 <= ((-30.787 + 1.1358 * X_N2Ha + 0.27297 * X_N2Ha**2)
                                                 * np.tanh(5.7409 * X_N2Ha) - 31.093))
    elif criterion == 'KA03':
        mask = (X_N2Ha <= (0.61 / (-1.7 - 1.3) + 0.05)) & \
               (OIII5006 <= (0.61 / (X_N2Ha - 0.05) + 1.3))
    elif criterion == 'KE01':
        mask = (X_N2Ha <= (0.61 / (-1.7 - 1.19) + 0.47)) & \
               (OIII5006 <= (0.61 / (X_N2Ha - 0.47) + 1.19))
    elif criterion == 'KE6A':
        mask = (X_N2Ha <= (0.61 / (-1.7 - 1.19) + 0.47)) & \
               (OIII5006 <= (0.61 / (X_N2Ha - 0.47) + 1.19)) & \
               (EW >= 6.)
    elif criterion == 'CF11':
        mask = (EW >= 3.) & (X_N2Ha <= -0.4)
    else:
        raise ValueError(f"Criterion '{criterion}' not recognized.")

    # --- Apply final mask ---
    x = np.array(r[mask])
    y = np.array(OH[mask])
    yerr = np.array(OH_err[mask])
    
    # --- Save CSV (optional) ---
    results = pd.DataFrame({
        'r': x,
        'OH': y,
        'eOH': yerr
    })
    
    if save_table:
    
        if calibrator == 1:
            calib = 'PP04_O3N2'
        elif calibrator == 2:
            calib = 'PP04_N2'
        elif calibrator == 3:
            calib = 'M13_O3N2'
        elif calibrator == 4:
            calib = 'M13_N2'
        elif calibrator == 5:
            calib = 'D16'
        else:
            raise ValueError("Invalid calibrator. Use 1=O3N2_PP04, 2=N2_PP04, 3=O3N2_M13, 4=N2_M13, 5=D16.")
            
        os.makedirs("tables_criterions", exist_ok=True)
    
        results.to_csv("tables_criterions/"+str(name)+"_"+str(criterion)+"_"+str(calib)+".csv", index=False)

    return x, y, yerr

