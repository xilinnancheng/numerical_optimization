import math
import random
import numpy as np
from matplotlib import pyplot as plt


def Rosenbrock(x):
    # Rosenbrock function
    f_val = 0
    half_x_size = int(x.size / 2)
    for i in range(half_x_size):
        f_val += 100 * pow(pow(x[2 * i], 2) - x[2*i+1],
                           2) + pow(x[2 * i] - 1, 2)
    return f_val


def Grad(x):
    # Rosenbrock gradient
    grad = []
    half_x_size = int(x.size / 2)
    for i in range(half_x_size):
        grad.append(400 * (pow(x[2 * i], 2) - x[2 * i + 1])
                    * x[2 * i] + 2*(x[2 * i] - 1))
        grad.append(200 * (x[2 * i + 1] - pow(x[2 * i], 2)))

    return np.array(grad)


if __name__ == "__main__":
    # Funtion dimention
    N = 2
    # Control parameter
    c = 0.7
    start_value = []
    for i in range(N):
        start_value.append(2.0 * 1)

    x = np.array(start_value)
    grad = Grad(x)
    x_list = np.array([x])
    iteration = 0

    # Main process
    while np.linalg.norm(grad) > 1e-5:
        tao = 1
        while Rosenbrock(x + tao * (- grad)) > Rosenbrock(x) + c * tao * np.dot(- grad, grad):
            tao = tao * 0.5
        x = x + tao * (- grad)
        grad = Grad(x)
        x_list = np.concatenate((x_list, [x]))
        iteration += 1
    print("iteration num:", iteration)
    print("optimal value x*=", x)

    # Plot Rosenbrock function
    if N == 2:
        plt.figure(1)
        x1 = np.arange(-2, 2, 0.05)
        x2 = np.arange(-1, 3, 0.05)
        x1, x2 = np.meshgrid(x1, x2)
        b = 100
        def f(x, y): return (x-1)**2 + b*(y-x**2)**2
        y = f(x1, x2)
        plt.contour(x1, x2, y, 200)
        plt.plot(x_list[:, 0], x_list[:, 1], "*", color='r')
        plt.title("Rosenbrock function")
        plt.show()
