#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.integrate import odeint
from pylab import * # for plotting commands
import numpy as np
import math
import random

print "3. Практическая часть"
print "3.1. Исследование для первого воздействия"
print "3.1.1. Построение имитатора объекта и ИИС"


# Установим объем выборки, предъявляемой к измерению
N = 1000
print "Объем выборки, предъявляемой к измерению: ", N

# Исходные данные
teta = [0.3, 0.2, 10.0]
print "Исходные данные: ", teta

# Входное воздействие
def x(t):
	if t >= 0:
		return 1
	else:
		return 0
print "Входное воздействие: "
print "\t1 if t >= 0"
print "\t0 else"

# Система дифференциальных уравнений, где Y[0] = y(t), Y[1] = y'(t)
def deriv(Y, t):
	return array([ Y[1], - teta[1] * Y[0] - teta[0] * Y[1] + teta[2] * x(t)])

print "Система дифференциальных уравнений, где Y[0] = y(t), Y[1] = y'(t): "
print "\tY[1]"
print "\t- teta[1] * Y[0] - teta[0] * Y[1] + teta[2] * x(t)"

# вектор начальных условий
yinit = array([1, 1])
print "Вектор начальных условий: ", yinit

# Временные границы эксперимента
tn = 0
tk = 100
print "Временные границы эксперимента: ", tn, tk
time = linspace(0.0, 100.0, N)

# Решаем получившееся ДУ
print "Решаем получившееся ДУ..."
y = odeint(deriv, yinit, time)
y1 = y[:, 0]
y2 = y[:, 1]
print "y1:\t\t\ty2:"
for i in xrange(5):
	print y1[i], "\t\t", y2[i]
print "...\t\t\t..."
for i in xrange(N - 5, N):
	print y1[i], "\t\t", y2[i]
print

# Рисуем график зависимости результатов моделирования от времени: у1 - выход модели, у2 - его производная
print "Рисуем график зависимости результатов моделирования от времени: у1 - выход модели, у2 - его производная (1.1.png)"
fig = figure()
subplot(111)
plot(time, y1, label="y1")
plot(time, y2, 'r--', label="y2")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
		ncol=2, mode="expand", borderaxespad=0.)
xlabel('t')
savefig("1.1.png")

# Построим имитатора объекта
print "Построим имитатора объекта"
# Предположим, что ошибки распределены по нормальному закону
n = np.random.normal(0, math.sqrt(0.01), N)

# Вектор измерений
Z = [y1[i] + n[i] for i in range(N)]

# Построим график имитатора объекта при наличии случайных возмущений на интервале 0..100 с delta = 0.1 с
print "Построим график имитатора объекта при наличии случайных возмущений на интервале 0..100 с delta = 0.1 с (1.2.png)"
fig.clear()
plot(time, Z)
xlabel('t')
ylabel('Z')
savefig("1.2.png")

# Построим график имитатора объекта при наличии случайных возмущений на интервале 0..2 с delta = 0.1 с
print "Построим график имитатора объекта при наличии случайных возмущений на интервале 0..100 с delta = 0.1 с (1.3.png)"
fig.clear()
time = linspace(0.0, 2.0, N / 50)
subplot(111)
plot(time, y1[:N / 50], label="y1")
plot(time, Z[:N / 50], 'r--', label="Z")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
		ncol=2, mode="expand", borderaxespad=0.)
xlabel('t')
savefig("1.3.png")

# Этап построение иммитатора объекта заканчивается на получении вектора измерений
print "Этап построение иммитатора объекта заканчивается на получении вектора измерений"
print "y1:\t\t\tn:\t\t\t\tZ:"
for i in xrange(5):
	print y1[i], "\t\t", n[i], "\t\t", Z[i]
print "...\t\t\t...\t\t\t\t..."
for i in xrange(N - 5, N):
	print y1[i], "\t\t", n[i], "\t\t", Z[i]
print

# 3.1.2. Проверка условий наблюдаемости и идентифицируемости объекта
print "3.1.2. Проверка условий наблюдаемости и идентифицируемости объекта\n"
# Проверим условие наблюдаемости:
print "Проверим условие наблюдаемости:"
# Выпишем матрицу А и матрицу измерений С
print "Матрица А:"
A = np.array([[0.0, 1.0], [-teta[1], -teta[0]]])
print A[0]
print A[1]
print "Матрица измерений C:"
C = np.array([1, 0])[np.newaxis]
print C[0]
# Сформируем матрицу наблюдаемости Калмана
print "Сформируем матрицу наблюдаемости Калмана:"
S = np.column_stack((C.T, A.T.dot(C.T)))
print S[0]
print S[1]
# Если ее ранг равен 2, то условие наблюдаемости выполняется
rank = np.linalg.matrix_rank(S)
if rank == 2:
	print "Ранг равен 2, отсюда следует, что объект является наблюдаемым"
else:
	print "Ранг равен ", rank, "отсюда следует, что объект не является наблюдаемым"

# Проверим условие идентифицируемости:
print "Проверим условие идентифицируемости:"

# Выпишем матрицу А и вектор начальных условий:
print "Матрица А:"
print A[0]
print A[1]
print "Вектор начальных условий:"
y0 = np.array([1, 1])
print y0

# Сформируем матрицу наблюдаемости Калмана
print "Сформируем матрицу наблюдаемости Калмана:"
S = np.column_stack((y0.T, A.T.dot(y0.T)))
print S[0]
print S[1]

rank = np.linalg.matrix_rank(S)
if rank == 2:
	print "Ранг равен 2, отсюда следует, что объект является идентифицируемым"
else:
	print "Ранг равен ", rank, "отсюда следует, что объект не является идентифицируемым"
print

# 3.1.3. Вычисление матрицы измерений С
print "3.1.3. Вычисление матрицы измерений С"
print
print "Составим и решим методом Рунге-Кутта систему уравнений чувствительности:"

print "Зададим матрицу частных производных:"

def deriv2(Y, t):
	return array([  Y[1],
					- teta[0] * Y[1] - teta[1] * Y[0] + teta[2] * x(t),
					Y[4],
					Y[5],
					- teta[1] * Y[2] - teta[0] * Y[4] - Y[1],
					- teta[1] * Y[3] - teta[0] * Y[5] - Y[0]])

print "\tY[1]"
print "\t- teta[0] * Y[1] - teta[1] * Y[0] + teta[2] * x(t)"
print "\tY[4]"
print "\tY[5]"
print "\t- teta[1] * Y[2] - teta[0] * Y[4] - Y[1]"
print "\t- teta[1] * Y[3] - teta[0] * Y[5] - Y[0]"

# вектор начальных условий
print "Зададим вектор начальных условий:"
U0 = array([1, 1, 0, 0, 0, 0])
print U0

# Временные границы эксперимента
time = linspace(0.0, 100.0, N)

# Решаем получившееся ДУ
print "Решение системы ДУ проводится с помощью функции odeint() посредством которой, выводится матрица U, столбцы которой содержат значения решений и производные."
U = odeint(deriv2, U0, time)
print "U:"
print U
# Присвоим интересующим нас параметрам значения матрицы U:
y1a = U[:, 0]
C1 = np.array([U[:, 2]])
C2 = np.array([U[:, 3]])

# Получим матрицу измерений:
print "Получим матрицу измерений из параметров матрицы U:"
C = np.column_stack((C1.T, C2.T))
print "C:"
print C

# Построим график чувствительности dy/d*teta1
print "Построим график чувствительности dy/d*teta1 (1.4.png)"
fig.clear()
time = linspace(0.0, 100.0, N)
plot(time, C[:, 0][:])
xlabel('t')
ylabel('C1')
savefig("1.4.png")

# Построим график чувствительности dy/d*teta2
print "Построим график чувствительности dy/d*teta2 (1.5.png)"
fig.clear()
time = linspace(0.0, 100.0, N)
plot(time, C[:, 1][:])
xlabel('t')
ylabel('C2')
savefig("1.5.png")

# 3.1.4. Реализация МНК-алгоритма
print "3.1.4. Реализация МНК-алгоритма"

# Определим требуемые оценки параметров
print "Определим требуемые оценки параметров:"
# teta_mnk = (C.T * C)^-1 * C.T * Z
print "Q = C.T * C:"
Q = np.matrix(C.T.dot(C)) # C.T * C
print Q
print "A = Q^-1"
A = Q.I # (C.T * C)^-1 
print A
# [  3.11617423e-06  -1.57025125e-08]
# [ -1.57025125e-08   1.60457589e-08]
print "b = C.T * Z:"
b = C.T.dot(Z)
print b
print "teta_mnk = A * b"
teta_mnk = A.dot(b)
teta_mnk = np.array([abs(teta_mnk[:, i]) for i in range(2)])
print teta_mnk

# 3.1.5. Построение эллипса рассеяния

# det(lambda*E - A)

# [lambda - 3.11617423e-06           1.57025125e-08]
# [         1.57025125e-08  lambda - 1.60457589e-08]

# (lamda - 3.11617423e-06)*(lambda - 1.60457589e-08) - 1.57025125e-08 * 1.57025125e-08 =
# = lambda^2 - 1.60457589e-08*lambda - 3.11617423e-06*lambda + 5.0001380384973146e-14 - 2.4656889881265625e-16=
# = lambda^2 - 3.1322199889e-06*lambda + 4.975481148616049e-14

#Подсчет и вывод корней уравнения    
def equation_solve(a,b,c):
    D=b**2-4*a*c
    if a:
        if D>0:
            x1=(-b+D**0.5)/(2*a)
            x2=(-b-D**0.5)/(2*a)
            print "Корни уравнения:\n","x1 =",x1,"\nx2 =",x2
            return [x1, x2]
        if D==0:
            x1=(-b)/(2*a)
            print "Корень уравнения:\n","x1 = x2 =",x1
        if D<0:
            print "Корни уравнения:"
            print 'x1 = '-b/2*a+math.sqrt(D)/2*abs(a)
            print 'x2 = '-b/2*a-math.sqrt(D)/2*abs(a)
    elif b:
        x1=-c/b
        print'Корень уравнения:\n','x =',x1
    elif c:
        print'Уравнение неверно'
    else:
        print'Уравнение верно'

# Находим длины полуосей
lambda1, lambda2 = equation_solve(1, - 3.1322199889e-06, 4.975481148616049e-14)
a1 = math.sqrt(lambda1) # 0.00176529141023
a2 = math.sqrt(lambda2) # 0.000126357531956

# Находим напрвление полуосей

# Для lambda1 = 3.11625376302e-06:
# [3.11625376302e-06 - 3.11617423e-06                      1.57025125e-08] [f1]   [0]
# [                    1.57025125e-08  3.11625376302e-06 - 1.60457589e-08] [f2] = [0]

# [7.953302000012681e-11     1.57025125e-08] [f1]   [0]
# [       1.57025125e-08  3.10020800412e-06] [f2] = [0]

# f1 = 1
# f2 = - 197.43387714907547

# Для lambda2 = 1.59662258821e-08:
# [1.59662258821e-08 - 3.11617423e-06                      1.57025125e-08] [f1]   [0]
# [                    1.57025125e-08  1.59662258821e-08 - 1.60457589e-08] [f2] = [0]

# [-3.1002080041179e-06          1.57025125e-08] [f1]   [0]
# [      1.57025125e-08  -7.953301790000126e-11] [f2] = [0]

# f1 = 1
# f2 = 0.005064986759284622

# Строим эллипс рассеяния
npts = 250
theta = np.arange(npts)*2.0*math.pi/(npts-1)
angle = np.arctan(- 197.43387714907547)
x = 0 + 0.00176529141023*np.cos(theta)*np.cos(angle) - 0.000126357531956*np.sin(theta)*np.sin(angle)
y = 0 + 0.00176529141023*np.cos(theta)*np.sin(angle) - 0.000126357531956*np.sin(theta)*np.cos(angle)
clf()
plot(x,y,color="r")
axis([-0.002, 0.002, -0.002, 0.002])
grid(True)
savefig("1.ellipse.png")
