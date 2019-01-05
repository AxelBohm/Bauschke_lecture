import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


class Line():

    def __init__(self, point, vector):
        self.point = point
        self.vec = vector
        self.unit_vec = vector / la.norm(vector)

    def projection(self, x):
        return np.inner(self.unit_vec, x-point) * self.unit_vec + point

    def resolvent(self, x):
        return 2*self.projection(x)-x


def resolvent_orthant(x):
    return 2*projection_orthant(x)-x


def projection_orthant(x):
    y = x.copy()
    for i in range(len(x)):
        y[i] = max(x[i], 0)

    return y


def alt_projec(x_0, line, iterations):
    x = x_0.copy()
    iter_save = []

    for _ in range(iterations):
        x = projection_orthant(line.projection(x))
        iter_save.append(x)

    return (x, iter_save)


def dr_operator(x, line):
    return 0.5*(x + line.resolvent(resolvent_orthant(x)))


def dr_algo(x_0, line, iterations):
    x = x_0.copy()
    iter_save = []
    shadow_save = []
    for _ in range(iterations):
        x = dr_operator(x, line)
        # store iterates
        iter_save.append(x)
        shadow_save.append(projection_orthant(x))

    return (x, iter_save, shadow_save)


# params
point = np.array([0.1, 0])
vector = np.array([0.1, -0.1])
x_0 = np.array([4, -0.1])
# --------

line_trough_int = Line(point, vector)

(x, iter_save, shadow_save) = dr_algo(x_0, line_trough_int, 15)
(x_alt_proj, iter_save_alt_proj) = alt_projec(x_0, line_trough_int, 15)

# prepare figure
fig_trough_int = plt.figure()
x_plot_iter = [xi[0] for xi in iter_save]
y_plot_iter = [xi[1] for xi in iter_save]
x_plot_shadow = [xi[0] for xi in shadow_save]
y_plot_shadow = [xi[1] for xi in shadow_save]
x_plot_alt_proj = [xi[0] for xi in iter_save_alt_proj]
y_plot_alt_proj = [xi[1] for xi in iter_save_alt_proj]
xlim = [min([*x_plot_iter, *x_plot_shadow, 0, *abs(x_0)])-1,
        max(*x_plot_iter, *x_plot_shadow, *abs(x_0), 1)+1]
ylim = [min(*y_plot_iter, *y_plot_shadow, 0, *abs(x_0))-1,
        max(*y_plot_iter, *y_plot_shadow, *abs(x_0))+1]


# plot line
slope = vector[1]/vector[0]
shift = point[1] - point[0]*slope
x_line = np.array(xlim)
y_line = x_line*slope + shift
plt.plot(x_line, y_line)

# plot iterates
plt.scatter(x_plot_iter, y_plot_iter, label='iterates')

# plot shadow sequence
plt.scatter(x_plot_shadow, y_plot_shadow, label='shadow sequence')

# plot nonnegative orthant
ylim[1] = max(*y_plot_iter, *y_plot_shadow, *y_line)+1
ax = plt.gca()
ax.fill_between([0, xlim[1]], 0, ylim[1], alpha=0.5)

# plot starting point
plt.scatter(*x_0, color='red', label='starting point')

plt.grid()
plt.legend(loc='lower left')
plt.savefig('intersec_int.eps')

# --------------------------------------------
# Next plot
# params
point = np.array([-1, 0])
vector = np.array([0.1, -0.1])
x_0 = np.array([4, -0.1])
# --------

line_trough_int = Line(point, vector)

(x, iter_save, shadow_save) = dr_algo(x_0, line_trough_int, 20)

# prepare figure
fig_trough_int = plt.figure()
x_plot_iter = [xi[0] for xi in iter_save]
y_plot_iter = [xi[1] for xi in iter_save]
x_plot_shadow = [xi[0] for xi in shadow_save]
y_plot_shadow = [xi[1] for xi in shadow_save]
xlim = [min([*x_plot_iter, *x_plot_shadow, 0, *abs(x_0)])-1,
        max(*x_plot_iter, *x_plot_shadow, *abs(x_0), 1)+1]
ylim = [min(*y_plot_iter, *y_plot_shadow, 0, *abs(x_0))-1,
        max(*y_plot_iter, *y_plot_shadow, *abs(x_0))+1]


# plot line
slope = vector[1]/vector[0]
shift = point[1] - point[0]*slope
x_line = np.array(xlim)
y_line = x_line*slope + shift
plt.plot(x_line, y_line)

# plot iterates
plt.scatter(x_plot_iter, y_plot_iter, label='iterates')

# plot shadow sequence
plt.scatter(x_plot_shadow, y_plot_shadow, label='shadow sequence')

# plot nonnegative orthant
ylim[1] = max(*y_plot_iter, *y_plot_shadow, *y_line)+1
ax = plt.gca()
ax.fill_between([0, xlim[1]], 0, ylim[1], alpha=0.5)

# plot starting point
plt.scatter(*x_0, color='red', label='starting point')

plt.grid()
plt.legend(loc='lower left')
plt.savefig('no_solution.eps')
