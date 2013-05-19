#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.integrate import odeint
from pylab import * # for plotting commands
import numpy as np
import random


# Установим объем выборки, предъявляемой к измерению
N = 1000

# Исходные данные
teta = [0.3, 0.2, 10]

# Входное воздействие
def x(t):
	if t >= 0:
		return 1
	else:
		return 0

# Система дифференциальных уравнений, где Y[0] = y(t), Y[1] = y'(t)
def deriv(Y, t):
	return array([ Y[1], -teta[1] * Y[0] - teta[0] * Y[1] + teta[2] * x(t)])

# вектор начальных условий
yinit = array([1, 1])

# Временные границы эксперимента
time = linspace(0.0, 10.0, N)

# Решаем получившееся ДУ
y = odeint(deriv, yinit, time)
y1 = y[:, 0]
y2 = y[:, 1]

# Рисуем график зависимости результатов моделирования от времени: у1- выход модели, у2 - его производная
fig = figure()
subplot(111)
plot(time, y1, label="y1")
plot(time, y2, label="y2")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
		ncol=2, mode="expand", borderaxespad=0.)
xlabel('t')
savefig("f1.png")

# Построим имитатора объекта

# Предположим, что ошибки распределены по нормальному закону
n = np.random.normal(0, 0.01, N)

# Вектор измерений
Z = [y1[i] + n[i] for i in range(N)]

# Построим график имитатора объекта при наличии случайных возмущений на интервале 0..10 с delta = 0.01 с
fig.clear()
plot(time, Z)
xlabel('t')
ylabel('Z')
savefig("f2.png")

# Построим график имитатора объекта при наличии случайных возмущений на интервале 0..1 с delta = 0.01 с
fig.clear()
time = linspace(0.0, 1.0, 100)
subplot(111)
plot(time, y1[:100], label="y1")
plot(time, Z[:100], label="Z")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
		ncol=2, mode="expand", borderaxespad=0.)
xlabel('t')
savefig("f3.png")
