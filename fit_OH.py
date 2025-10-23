import distance
import abundance
import criteria
import models
import plot

def fit_final(name, HIIREGID, ra, ra0, dec, dec0, pa, ba, d, re, EWHa, Hb4861, eHb4861, Ha6562, eHa6562, OIII5006, eOIII5006, NII6583, eNII6583, SII6716, eSII6716, SII6730, eSII6730, calibrator, criterion, save_table, save_graph, show_graph):

    x = distance.distances(ra, ra0, dec, dec0, pa, ba, d, re)
    
    y, ey, Ha6562_cor, OIII5006_cor, NII6583_cor = abundance.abundance(name, x, HIIREGID, EWHa, Hb4861, eHb4861, Ha6562, eHa6562, OIII5006, eOIII5006, NII6583, eNII6583, SII6716, eSII6716, SII6730, eSII6730, calibrator)
    
    r, oh, eoh = criteria.points(name, criterion, x, y, ey, EWHa, Ha6562_cor, OIII5006_cor, NII6583_cor, calibrator, save_table)
    
    results_dict = models.fit_models(r, oh, eoh)

    if results_dict is not None:
        output = plot.plot_model(results_dict, name, criterion, calibrator, save_graph, show_graph)
        print(f"{name} Completed!")  # this will be displayed
        return output
    else:
        print(f"Insufficient data for fitting the galaxy {name}.")
        return None
