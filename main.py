import matplotlib.pyplot as plt
import numpy as np
from random import randint

n = [10, 20]  # Число сплайнов
m = [2, 3, 4]
flag = [True, False]
begin, end, k = -1, 2, 0.001  # Интервал


# Полином 3 степени
def poly(a, b, c, d, x, x0):
    return a + b * (x - x0) + 1 / 2 * c * (x - x0) ** 2 + 1 / 6 * d * (x - x0) ** 3


# Рисование сетки и осей
def drawAxesGrid(color='red', width=0.8):
    plt.axhline(y=0, color=color, linewidth=width)
    plt.axvline(x=0, color=color, linewidth=width)
    plt.grid()


# Генерация системы для среднеквадратичного приближения
def MakeSystem(xyTable, basis):
    matrix = [[0] * (basis) for _ in range(basis)]
    right = [0] * (basis)
    for i in range(basis):
        for j in range(basis):
            sumA, sumB = 0, 0
            for k in range(len(xyTable)):
                sumA += xyTable[k][0] ** (i + j)  # * xyTable[k][0] ** j
                sumB += xyTable[k][1] * xyTable[k][0] ** i

            matrix[i][j] = sumA
            right[i] = sumB

    return matrix, right


# Полином для среднеквадратичного приближения
def polyC(c, x):
    sum = 0
    for i in range(len(c)):
        sum += c[i] * x ** i
    return sum


# Генерация шума
def rand(a, b, k):
    p = 0.1
    l = int((b - a) / k * p)
    for _ in range(5):
        n = (randint(0, l) - 0.5 * l) * k
        return n


# Генерация x на интервале [a, b]
def randX(a, b):
    x = (randint(0, int(b - a) * 1000) - 1000) / 1000
    return x


# Сглаживание результатов эксперимента
def sliding(y):
    for i in range(len(y) - 1):
        y[i] = (y[i] + y[i - 1] + y[i + 1]) / 3
    y[0] = (5 * y[0] + 2 * y[1] - y[2]) / 6
    y[- 1] = (5 * y[-1] + 2 * y[-2] - y[-3]) / 6
    return y


# Исходная функция
def func(x):
    return (5 * np.sin(5 * x) + x) / (np.cos(x) + 2 * x - 5)


def xSearch(begin, end, n, flag=True):
    if flag:
        # Вычисление x
        x = []
        for i in range(n + 1):
            xTmp = begin + (end - begin) * i / n
            x.append(xTmp)
    else:
        x = [begin]
        for i in range(n - 1):
            xTmp = randX(begin, end)
            x.append(xTmp)
        x.append(end)
        x.sort()
    return x


def ySearch(x):
    y = []

    for j in range(len(x)):
        y.append(func(x[j]) + rand(begin, end, k))
    y = sliding(y)
    return y


# ______________________________________
# ПОСТРОЕНИЕ СРЕДНЕКВАДРАТИЧНОГО ПРИБЛИЖЕНИЯ
# ______________________________________
def squareSearch(n, m, begin, end):
    x = xSearch(begin, end, n)
    y = ySearch(x)
    xy = []
    for k in range(len(x)):
        xy.append([x[k], y[k]])
    matrix, right = MakeSystem(xy, m)
    cPoly = np.linalg.solve(matrix, right)
    return x, y, cPoly


# ______________________________________
# ПОСТРОЕНИЕ СПЛАЙНОВ
# ______________________________________
def splineSearch(begin, end, n, flag):
    x = xSearch(begin, end, n, flag)

    # Вычисление y
    y = ySearch(x)

    # Вычисление  коэффициента a
    a = []
    for j in range(len(y)):
        a.append(y[j])

    h = [0]  # Вычисление h (шаг интерполирования)
    for j in range(1, len(x)):
        h.append(x[j] - x[j - 1])

    # Составление системы уравнений для вычисления коэффициента c
    A = [[0] * n for i in range(n)]
    A[0][0] = 2 * (h[1] + h[2])
    A[0][1] = h[2]
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            if i == j:
                A[i][j - 1] = h[i + 1]
                A[i][j] = 2 * (h[i + 1] + h[i + 2])
                A[i][j + 1] = h[i + 2]
    A.pop(-1)
    for j in range(len(A)):
        A[j].pop(-1)

    rightMatrix = [0]
    for j in range(1, n):
        rTmp = 6 * ((y[j + 1] - y[j]) / h[j + 1] - (y[j] - y[j - 1]) / h[j])
        rightMatrix.append(rTmp)
    rightMatrix.pop(0)

    #  Вычисление коэффициента c
    c = [0]
    cTmp = np.linalg.solve(A, rightMatrix)
    for j in range(len(cTmp)):
        c.append(cTmp[j])
    c.append(0)

    #  Вычисление коэффициента d
    d = [0]
    for j in range(1, len(c)):
        dTmp = (c[j] - c[j - 1]) / h[j]
        d.append(dTmp)

    #  Вычисление коэффициента b
    b = [0]
    for j in range(1, len(c)):
        bTmp = 0.5 * h[j] * c[j] - h[j] ** 2 * d[j] / 6 + (y[j] - y[j - 1]) / h[j]
        b.append(bTmp)

    return x, y, a, b, c, d


# ______________________________________
# РИСОВАНИЕ СПЛАЙНОВ
# ______________________________________
def splineDraw(a, b, c, d, x, n):
    for j in range(1, len(a)):

        xs = np.arange(x[j - 1], x[j] + k, k)
        if j == 1:
            plt.plot(xs, poly(a[j], b[j], c[j], d[j], xs, x[j]), color='red', label='Сплайны (n=' + str(n) + ')')
        else:
            plt.plot(xs, poly(a[j], b[j], c[j], d[j], xs, x[j]), color='red')


for i in range(len(n)):
    for j in range(len(flag)):
        x, y, a, b, c, d = splineSearch(begin, end, n[i], flag[j])

        plt.figure(figsize=(10, 10))
        plt.title(
            "Графики исходной функции и кубических сплайнов сплайнов \n Число узлов сетки (n = " + str(n[i]) + ")")

        plt.scatter(x, y, color='green', label='Экспериментальные данные')
        splineDraw(a, b, c, d, x, n[i])
        xInterval = np.arange(begin, end + k, k)
        plt.plot(xInterval, func(xInterval), color='blue', label='Исходная функция')
        drawAxesGrid('black')
        plt.legend()

        if flag[j]:
            string = 'regular'
            # plt.savefig('spline_' + str(n[i]) + '_ regular.png', bbox_inches='tight')
        else:
            string = 'irregular'
            # plt.savefig('spline_' + str(n[i]) + '_ irregular.png', bbox_inches='tight')
        plt.savefig('image/spline/spline_' + str(n[i]) + '_ ' + string + '.png', bbox_inches='tight')

# _____________________________
# РИСОВАНИЕ СРЕДНЕКВАДРАТИЧНОГО ПРИБЛИЖЕНИЯ
# _____________________________

for i in range(len(n)):
    m.append(n[i] - 1)
    for j in range(len(m)):
        x, y, cPoly = squareSearch(n[i], m[j], begin, end)
        xs = np.arange(x[0], x[-1], 0.01)

        plt.figure(figsize=(10, 10))
        plt.title(
            "Графики исходной функции и среднеквадратичного приближения "
            "\n Число узлов сетки (n = " + str(n[i]) + ")"
                                                       "\n Степень полинома (m = " + str(m[j]) + ")")

        xInterval = np.arange(begin, end + k, k)
        plt.plot(xInterval, func(xInterval), color='blue', label='Исходная функция')

        plt.plot(xs, polyC(cPoly, xs), color='red', label='Среднеквадратичное приближение (m=' + str(m[j]) + ')')

        plt.scatter(x, y, color='green', label='Экспериментальные данные')

        drawAxesGrid('black')
        plt.legend()
        plt.savefig('image/square/square_' + str(n[i]) + '_' + str(m[j]) + '.png', bbox_inches='tight')
    m.pop()
