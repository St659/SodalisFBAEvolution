import matplotlib.pyplot as plt
import glob
from Precap_Gloss_Heatmap import natural_sort_key
import numpy as np
import os
import csv
import matplotlib.cm as cm

results_directory = './iRHx_blood'
fig, growth_plot = plt.subplots()
fig2, hist_plot = plt.subplots()
results_files = glob.glob1(results_directory,'*.npz')
print(results_files)

# colors = ['xkcd:hot pink','xkcd:pink','xkcd:dark pink','xkcd:light pink']

colours = cm.winter(np.linspace(0,1,len(results_files)))
with open ('iRHx_blood_reactions.csv','w') as csv_a:
    field_names = ['growth', 'num_reactions', 'reaction_names']
    writer_a = csv.DictWriter(csv_a, fieldnames=field_names)
    # writer_b = csv.DictWriter(csv_b, fieldnames=field_names)

    for file, colour in zip(sorted(results_files, key=natural_sort_key), colours):

        print(file)
        saved_results = np.load(os.path.join(results_directory,file))
        # hist_plot.hist(saved_results['growth'],bins=1000)
        for growth, reaction_num, names in zip(saved_results['growth'], saved_results['used_reactions'],
                                               saved_results['used_reaction_names']):
            if growth < 0.42 and growth > 0.18 and reaction_num < 130:
                print('Growth: ' + str(growth) + ' Reactions: ' + str(reaction_num))
                print(names)
                writer_a.writerow({'growth': growth, 'num_reactions': reaction_num, 'reaction_names': names})
        #     if growth < 0.3 and growth > 0.2 and reaction_num < 110:
        #         print('Growth: ' + str(growth) + ' Reactions: ' + str(reaction_num))
        #         print(names)
        #         writer_b.writerow({'growth': growth, 'num_reactions': reaction_num, 'reaction_names': names})
        print(len(saved_results['used_reactions']))
        growth_plot.plot(saved_results['used_reactions']+229,saved_results['growth'],'o',color=colour)
growth_plot.set_xlabel('Reactions')
growth_plot.set_ylabel('Biomass Output')
# fig.savefig('./Glucose.png',dpi=1200)
growth_plot.set_title('iRHx_blood')
fig.savefig('./paper_figures/blood.pdf',format='pdf')
# plt.show()