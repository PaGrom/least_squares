#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.integrate import odeint
from pylab import * # for plotting commands


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

# Рисуем график зависимости результатов моделирования от времени: у1- выход модели, у2 - его производная
figure()
subplot(111)
plot(time, y[:, 0], label="y1")
plot(time, y[:, 1], label="y2")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
		ncol=2, mode="expand", borderaxespad=0.)
xlabel('t')
savefig("f1.png")
show()