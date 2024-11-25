# -*- coding: utf-8 -*-
import sys
import pandas as pd
import matplotlib.pyplot as plt

datafile = sys.argv[1]


data = pd.read_csv(datafile)

max_days = data['Time'].max()

if max_days >= 3650/2: 
    time_scale = 'Years'
    data['Time_Scaled'] = data['Time'] / 365
elif max_days >= 365: 
    time_scale = 'Months'
    data['Time_Scaled'] = data['Time'] / 30
elif max_days >= 7:  
    time_scale = 'Weeks'
    data['Time_Scaled'] = data['Time'] / 7
else:
    time_scale = 'Days'
    data['Time_Scaled'] = data['Time']


plt.figure(figsize=(18, 6))


plt.subplot(1, 3, 1)
plt.plot(data['Time_Scaled'], data['Prey'], label='Prey', color='blue')
plt.xlabel(time_scale)
plt.ylabel('Population')
plt.title(f'{datafile} Prey Population Over Time')
plt.legend()
plt.grid(True)

plt.subplot(1, 3, 2)
plt.plot(data['Time_Scaled'], data['Predator'], label='Predators', color='red')
plt.xlabel(time_scale)
plt.ylabel('Population')
plt.title(f'{datafile} Predator Population Over Time')
plt.legend()
plt.grid(True)

plt.subplot(1, 3, 3)
plt.scatter(data['Prey'], data['Predator'], label='Prey vs Predators', color='green')
plt.xlabel('Prey Population')
plt.ylabel('Predator Population')
plt.title(f'{datafile} Prey vs Predator Population')
plt.legend()
plt.grid(True)

# Финальная настройка и отображение
plt.tight_layout()
plt.show()
