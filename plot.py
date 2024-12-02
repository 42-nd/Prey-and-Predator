# -*- coding: utf-8 -*-
import sys
import pandas as pd
import matplotlib.pyplot as plt

datafile = sys.argv[1]

data = pd.read_csv(datafile)

plt.figure(figsize=(18, 6))


plt.subplot(1, 3, 1)
plt.plot(data['Time'], data['Prey'], label='Prey', color='blue')
plt.xlabel('Days')
plt.ylabel('Population')
plt.title(f'{datafile} Prey Population Over Time')
plt.legend()
plt.grid(True)


plt.subplot(1, 3, 2)
plt.plot(data['Time'], data['Predator'], label='Predators', color='red')
plt.xlabel('Days')
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


plt.tight_layout()
plt.show()
