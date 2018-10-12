#!/usr/bin/env python2
# -*- coding: utf-8 -*-

''' see
http://benalexkeen.com/bar-charts-in-matplotlib/
'''

import numpy as np
import matplotlib.pyplot as plt

countries = ['USA', 'GB', 'China', 'Russia', 'Germany']
bronzes = np.array([38, 17, 26, 19, 15])
silvers = np.array([37, 23, 18, 18, 10])
golds = np.array([46, 27, 26, 19, 17])
ind = [x for x, _ in enumerate(countries)]

total = bronzes + silvers + golds
proportion_bronzes = np.true_divide(bronzes, total) * 100
proportion_silvers = np.true_divide(silvers, total) * 100
proportion_golds = np.true_divide(golds, total) * 100

plt.bar(ind, proportion_golds, width=0.8, label='golds', color='gold', bottom=proportion_bronzes+proportion_silvers)
plt.bar(ind, proportion_silvers, width=0.8, label='silvers', color='silver', bottom=proportion_bronzes)
plt.bar(ind, proportion_bronzes, width=0.8, label='bronzes', color='#CD853F')

plt.xticks(ind, countries)
plt.ylabel("Medals")
plt.xlabel("Countries")
plt.title("2012 Olympics Top Scorers' Medals by Proportion")
plt.legend(bbox_to_anchor=(1.04,0.8), loc="center left", borderaxespad=0)
plt.ylim=1.0

# rotate axis labels
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')

plt.show()