#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.integrate import odeint
from pylab import * # for plotting commands
import numpy as np
import math
import random

# 3. Практическая часть
# 3.1. Исследование для первого воздействия
# 3.1.1. Построение имитатора объекта и ИИС


# Установим объем выборки, предъявляемой к измерению
N = 1000

# Исходные данные
teta = [0.3, 0.2, 10.0]

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
time = linspace(0.0, 100.0, N)

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
n = np.random.normal(0, math.sqrt(0.01), N)

# Вектор измерений
Z = [y1[i] + n[i] for i in range(N)]

# Построим график имитатора объекта при наличии случайных возмущений на интервале 0..100 с delta = 0.01 с
fig.clear()
plot(time, Z)
xlabel('t')
ylabel('Z')
savefig("f2.png")

# Построим график имитатора объекта при наличии случайных возмущений на интервале 0..10 с delta = 0.01 с
fig.clear()
time = linspace(0.0, 10.0, 100)
subplot(111)
plot(time, y1[:100], label="y1")
plot(time, Z[:100], label="Z")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
		ncol=2, mode="expand", borderaxespad=0.)
xlabel('t')
savefig("f3.png")

# Этап построение иммитатора объекта заканчивается на получении вектора измерений
print "y1:\t\t\tn:\t\t\t\tZ:"
for i in xrange(15):
	print y1[i], "\t\t", n[i], "\t\t", Z[i]


# 3.1.2. Проверка условий наблюдаемости и идентифицируемости
# объекта

# Проверим условие наблюдаемости:

# Выпишем матрицу А и матрицу измерений С
A = np.array([[0.0, 1.0], [-teta[1], -teta[0]]])
C = np.array([1, 0])[np.newaxis]

# Сформируем матрицу наблюдаемости Калмана
S = np.column_stack((C.T, A.T.dot(C.T)))
# Если ее ранг равен 2, то условие наблюдаемости выполняется
rank = np.linalg.matrix_rank(S)
print rank
# rank = 2, отсюда следует, что объект является наблюдаемым

# Проверим условие идентифицируемости:

# Выпишем матрицу А и вектор начальных условий:
A = np.array([[0.0, 1.0], [-teta[1], -teta[0]]])
y0 = np.array([1, 1])

# Сформируем матрицу наблюдаемости Калмана
S = np.column_stack((y0.T, A.T.dot(y0.T)))

rank = np.linalg.matrix_rank(S)
print rank
# rank = 2, отсюда следует, что объект является идентифицируемым

# 3.1.3. Вычисление матрицы измерений С

# Зададим некие априорные значения параметров:

teta_a = [60, 0.005, 100]

def deriv2(Y, t):
	return array([ Y[1],
					- teta_a[0] * Y[1] - teta_a[1] * Y[0] + teta_a[2] * x(t),
					Y[4],
					Y[5],
					- teta_a[1] * Y[2] - teta_a[0] * Y[4] - Y[1],
					- teta_a[1] * Y[3] - teta_a[0] * Y[5] - Y[0]])

# вектор начальных условий
U0 = array([1, 1, 0, 0, 0, 0])

# Временные границы эксперимента
time = linspace(0.0, 100.0, N)

# Решаем получившееся ДУ
U = odeint(deriv2, U0, time)

# Присвоим интересующим нас параметрам значения матрицы U:
y1a = U[:, 0]
C1 = np.array([U[:, 2]])
C2 = np.array([U[:, 3]])

# Получим матрицу измерений:
C = np.column_stack((C1.T, C2.T))

# Построим график cравнение выхода модели у1 и вычисленного по априорным оценкам выхода модели на интервале 0..10 с
fig.clear()
time = linspace(0.0, 10.0, 100)
subplot(111)
plot(time, y1[:100], label="y1")
plot(time, y1a[:100], label="y1a")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
		ncol=2, mode="expand", borderaxespad=0.)
xlabel('t')
savefig("f4.png")

# Построим график чувствительности dy/d*teta1
fig.clear()
time = linspace(0.0, 10.0, 100)
plot(time, C[:, 0][:100])
xlabel('t')
ylabel('C1')
savefig("f5.png")

# Построим график чувствительности dy/d*teta2
fig.clear()
time = linspace(0.0, 10.0, 100)
plot(time, C[:, 1][:100])
xlabel('t')
ylabel('C2')
savefig("f6.png")

