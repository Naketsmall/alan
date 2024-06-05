import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

import converter


DS = pd.read_csv('/home/epsilon/Downloads/Telegram Desktop/t0_HSL+rates_corrected (2).dat', sep="\s+")

# вставили ручками

# из D_MW.COP
M = [0.001, 0.016, 0.017, 0.018, 0.032, 0.002, 0.028]

# из Therm.cop
Y_names = ['YH', 'YO', 'YOH', 'YH2O', 'YO2', 'YH2', 'YN2']
A = [211798519.7, 15188041.26, 1705091.865,-14030832.22, -315060.297, -4036848.415, -337602.0651]
B = [20786031.5, 1315517.835, 1621353.945, 1803438.791, 980518.5471, 13495527.88, 1058584.796]
C = [0.0007681104921, 3130.44411, 102031.2378, 257854.3523, 49804.08554, 882411.5031, 54318.53763]

eqs = [[0, 0, -2, 0, 1, 1, 0],
       [1, -1, -1, 0, 1, 0, 0],
       [-1, 0, 1, -1, 0, 1, 0],
       [-1, 1, -1, 0, 0, 1, 0],
       [0, -1, 2, -1, 0, 0, 0],
       [2, 0, 0, 0, 0, -1, 0],
       [1, 0, 1, -1, 0, 0, 0]]


conv = converter.Converter(DS, M, Y_names, A, B, C, eqs)

q1 = conv.calculate_q1()

#plt.plot(DS['x'], q1)
#plt.plot(DS['x'], DS['Q'])'''

#Q1 = conv.calculate_Q1()
#Q2 = conv.calculate_Q2()
Q4 = conv.calculate_Q4()
#Q3 = conv.calculate_Q3()

#dqQ1 = conv.calculate_dqQ(Q1)
#dqQ2 = conv.calculate_dqQ(Q2)
dqQ4 = conv.calculate_dqQ(Q4)
#dqQ3 = conv.calculate_dqQ(Q3)

#qQ1 = converter.integrate_var(DS['x'], dqQ1)
#qQ2 = converter.integrate_var(DS['x'], dqQ2)
qQ4 = converter.integrate_var(DS['x'], dqQ4)
print(qQ4)
#qQ3 = converter.integrate_var(DS['x'], dqQ3)

#plt.plot(DS['x'], q1, color='orange', label='q1')
#plt.plot(DS['x'], DS['Q'], color='red', label="DS['Q']")
#plt.plot(DS['x'], qQ1, color='blue', label='q(Q1)')
#plt.plot(DS['x'], qQ2, color='purple', label='q(Q2)')
#plt.plot(DS['x'], Q2[:,5:7], label='q(Q2)')
plt.plot(DS['x'], qQ4, label='q(Q2)')
#plt.plot(DS['x'], qQ3, color='green', label='q(Q3)')
#plt.plot(DS['x'], [dqQ2[i]-dqQ1[i] for i in range(len(DS))], color='red')
#plt.plot(DS['x'], Q2-Q1, color='blue')
#plt.plot(DS['x'], DS['p_total'], color='orange')
#plt.plot(DS['x'], DS['p'])
plt.legend(loc='best')
plt.grid()
plt.show()

