import matplotlib.pyplot as plt
import numpy as np
import os

def plot_model(results_dict, name, criterion, calibrator, save_graph, show_graph):
    
    x = results_dict['x']
    y = results_dict['y']
    ey = results_dict['ey']
    best_case = results_dict['best_case']

    fig, ax = plt.subplots()
    ax.errorbar(x, y, ey, fmt='o', color='black', elinewidth=0.5, 
                markerfacecolor='mediumslateblue', markeredgecolor='black', ms=5, capsize=3, zorder=1)
    
    ymin, ymax = 0.97*min(y), 1.03*max(y)
    
    if best_case == 1:
        a2 = results_dict['fit1'][0]
        ea2 = results_dict['fit1'][1]
        b0 = results_dict['fit1'][2]
        eb0 = results_dict['fit1'][3]
        p = np.linspace(min(x), max(x), 100)
        ax.plot(p, a2*p + b0, color='red', linewidth=2, zorder=3)
        
        b0 = b0
        eb0 = eb0
        a1 = ea1 = h1 = eh1 = 0.0
        a2 = a2
        ea2 = ea2
        h2 = eh2 = a3 = ea3 = 0.0
        
        print('h1 = {:.2f}'.format(h1), 'a1 = {:.2f}'.format(a1), 'b1 = {:.2f}'.format(b0))
    
    elif best_case == 2:
        fit2 = results_dict['fit2']
        fit2.plot_fit(color='red', linewidth=2, zorder=3)
        fit2.plot_breakpoints(linewidth=0.8, color='gray', linestyle='-.', zorder=2)
        
        interval0 = fit2.get_results()["estimates"]["breakpoint1"]["confidence_interval"][0]
        interval1 = fit2.get_results()["estimates"]["breakpoint1"]["confidence_interval"][1]
        rect1 = plt.Rectangle((interval0, ymin), (interval1-interval0), ymax-ymin, facecolor="azure", zorder=0.5)
        ax.add_patch(rect1)
        
        b0 = fit2.get_results()["estimates"]["const"]["estimate"]
        eb0 = fit2.get_results()["estimates"]["const"]["se"]
        a1 = fit2.get_results()["estimates"]["alpha1"]["estimate"]
        ea1 = fit2.get_results()["estimates"]["alpha1"]["se"]
        h1 = fit2.get_results()["estimates"]["breakpoint1"]["estimate"]
        eh1 = fit2.get_results()["estimates"]["breakpoint1"]["se"]
        a2 = fit2.get_results()["estimates"]["alpha2"]["estimate"]
        ea2 = fit2.get_results()["estimates"]["alpha2"]["se"]
        h2 = eh2 = a3 = ea3 = 0.0
        
        print('h1 = {:.2f}'.format(h1), 'a1 = {:.2f}'.format(a1), 'b1 = {:.2f}'.format(b0))
        
    else:
        fit3 = results_dict['fit3']
        fit3.plot_fit(color='red', linewidth=2, zorder=3)
        fit3.plot_breakpoints(linewidth=0.8, color='gray', linestyle='-.', zorder=2)    
        
        interval0 = fit3.get_results()["estimates"]["breakpoint1"]["confidence_interval"][0]
        interval1 = fit3.get_results()["estimates"]["breakpoint1"]["confidence_interval"][1]
        rect1 = plt.Rectangle((interval0, ymin), (interval1-interval0), ymax-ymin, facecolor="azure", zorder=0.5)
        ax.add_patch(rect1)
        
        interval00 = fit3.get_results()["estimates"]["breakpoint2"]["confidence_interval"][0]
        interval11 = fit3.get_results()["estimates"]["breakpoint2"]["confidence_interval"][1]
        rect2 = plt.Rectangle((interval00, ymin), (interval11-interval00), ymax-ymin, facecolor="azure", zorder=0.5)
        ax.add_patch(rect2)
        
        b0 = fit3.get_results()["estimates"]["const"]["estimate"]
        eb0 = fit3.get_results()["estimates"]["const"]["se"]
        a1 = fit3.get_results()["estimates"]["alpha1"]["estimate"]
        ea1 = fit3.get_results()["estimates"]["alpha1"]["se"]
        h1 = fit3.get_results()["estimates"]["breakpoint1"]["estimate"]
        eh1 = fit3.get_results()["estimates"]["breakpoint1"]["se"]
        a2 = fit3.get_results()["estimates"]["alpha2"]["estimate"]
        ea2 = fit3.get_results()["estimates"]["alpha2"]["se"]
        h2 = fit3.get_results()["estimates"]["breakpoint2"]["estimate"]
        eh2 = fit3.get_results()["estimates"]["breakpoint2"]["se"]
        a3 = fit3.get_results()["estimates"]["alpha3"]["estimate"]
        ea3 = fit3.get_results()["estimates"]["alpha3"]["se"]
        
        print('h1 = {:.2f}'.format(h1), 'a1 = {:.2f}'.format(a1), 'b1 = {:.2f}'.format(b0))
        
    ax.set_title(str(name))
    ax.set_xlabel("$r$/r$_e$", size=14)
    ax.set_ylabel("12+log(O/H)", size=14)
    ax.set_xlim([0, max(x)*1.05])
    ax.set_ylim([ymin, ymax])
    ax.minorticks_on()
    ax.tick_params(which='major', direction='in', length=4.0, width=0.7, colors='black', grid_color='gray', grid_alpha=0.9)
    ax.tick_params(which='minor', direction='in', length=2.0, width=0.5, colors='black', grid_color='gray', grid_alpha=0.9)
    
    if save_graph:
    
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
            
        os.makedirs("graphs", exist_ok=True)
    
        plt.savefig("graphs/"+str(name)+"_"+str(criterion)+"_"+str(calib)+".png", transparent=False, facecolor='w', edgecolor='w')
        
    # Show the plot or not
    if show_graph:
        plt.show(block=True)
    else:
        plt.close(fig)
        
    return {
        'galaxy': name,
        'b0': b0,
        'eb0': eb0,
        'a1': a1,
        'ea1': ea1,
        'h1': h1,
        'eh1': eh1,
        'a2': a2,
        'ea2': ea2,
        'h2': h2,
        'eh2': eh2,
        'a3': a3,
        'ea3': ea3
    }

