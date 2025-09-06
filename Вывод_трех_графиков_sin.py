import matplotlib.pyplot as plt
import numpy as np

# Загрузка данных
data_godunov2 = np.loadtxt('results_G2.txt', delimiter='\t')
data_godunov1 = np.loadtxt('results_G1.txt', delimiter='\t')

x = data_godunov2[:, 0]
q_godunov1 = data_godunov1[:, 1]
q_godunov2 = data_godunov2[:, 1]
q_exact =  data_godunov2[:, 2]

# Находим точку разрыва
grad = np.abs(np.gradient(q_godunov2))
discontinuity_idx = np.argmax(grad)
discontinuity_x = x[discontinuity_idx]

# Построение графиков
plt.figure(figsize=(12, 6))
plt.plot(x, q_godunov1, 'bo--', markersize=4, linewidth=1, label='Численное решение по схеме Годунова 1-го порядка')
plt.plot(x, q_godunov2, 'go--', markersize=4, linewidth=1, label='Численное решение по схеме Годунова 2-го порядка')
plt.plot(x, q_exact, 'r-', linewidth=2, label='Точное решение')

# Отмечаем разрыв на оси
plt.axvline(x=discontinuity_x, color='k', linestyle=':', linewidth=1)
plt.text(discontinuity_x, plt.ylim()[0], f' x={discontinuity_x:.2f}', 
         color='r', ha='left', va='bottom')


plt.xlabel('Position', fontsize=12)
plt.ylabel('q', fontsize=12)
plt.title('Сравнение схем Годунова 1-го и 2-го порядков с точным решением (t = 1.0, q = sin(2*pi*x/1.0))', fontsize=14)
plt.legend()
plt.grid(True)
plt.savefig('comparison.png', dpi=300)
plt.show()
