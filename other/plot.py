import matplotlib.pyplot as plt

def plot(x,y, title, x_axis_label = 'x', y_axis_label = 'y', line_style = 'b-', xscale = 'log', yslace = 'log'):

    plt.plot(x, y, line_style, linewidth=3)
    plt.xscale(xscale)
    plt.yscale(yslace)
    plt.rcParams.update({'font.size': 16})
    plt.ylabel(y_axis_label)
    plt.xlabel(x_axis_label)
    plt.title(title)
    plt.xlim(1.,6.5)
    plt.show()

def double_plot(x_1, y_1, x_2, y_2, title,label_switch = 0, y_1_label = None, y_2_label = None,
                x_axis_label = 'x', y_axis_label = 'y', line_style_1 = 'b-', line_style_2 = 'y-', xslace = 'log', yslace = 'log'):

    plt.xscale(xslace)
    plt.yscale(yslace)
    plt.rcParams.update({'font.size': 16})
    plt.ylabel(y_axis_label)
    plt.xlabel(x_axis_label)
    plt.title(title)

    if (label_switch):

        plt.plot(x_1 , y_1, line_style_1, linewidth=2, label=y_1_label)
        plt.plot(x_2 , y_2, line_style_2, linewidth=2, label=y_2_label)
        plt.legend(loc=4,prop={'size':16})

    else:

        plt.plot(x_1 , y_1, line_style_1, linewidth=3)
        plt.plot(x_2 , y_2, line_style_2, linewidth=3)

    plt.show()


def multi_plot(X, Y, title, labels, label_switch = 0,
               x_axis_label = 'x', y_axis_label = 'y', line_styles = 'b-', xslace = 'log', yslace = 'log'):

    plt.xscale(xslace)
    plt.yscale(yslace)
    plt.rcParams.update({'font.size': 16})
    plt.ylabel(y_axis_label)
    plt.xlabel(x_axis_label)
    plt.title(title)

    if (label_switch):

        for i in range(0,len(Y[0,:])):

            plt.plot(X[:,i] , Y[:,i], line_styles[i], linewidth=2, label=labels[i])

        plt.legend(loc=3,prop={'size':8})

    else:

        for i in range(0,len(Y[0,:])):

            plt.plot(X[:,i] , Y[:,i], line_styles[i], linewidth=3)
    plt.show()

