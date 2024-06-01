import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import plotly.express as px
import re
import numpy as np
import os.path

def plot_projectile(trajectory_path: str, title: str, save_name: str = None):

    earth = pd.read_csv('satellite/trajectory/earth.csv')
    points = pd.read_csv(trajectory_path)

    points = points.iloc[::1000]
    R = 6378.137

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(earth.x.values, earth.y.values, earth.z.values, 'green')
    ax.plot(points.x.values, points.y.values, points.z.values, color='red')
    ax.set(xlabel='x-coordinate', ylabel='y-coordinate', zlabel='z-coordinate', title=title)
    ax.grid()
    fig.tight_layout()

    if isinstance(save_name, str):
        fig.savefig(fname=save_name, bbox_inches='tight')

    name, ext = os.path.splitext(save_name)
    coords_name = [(points.x, points.y), (points.x, points.z), (points.y, points.z)]
    earth_coords_name = [(earth.x, earth.y), (earth.x, earth.z), (earth.y, earth.z)]
    plots_2d = [('x', 'y'), ('x', 'z'), ('y', 'z')]
    for i in range(len(plots_2d)):
        fig = plt.figure()
        ax = fig.add_subplot()
        # ax.plot(*earth_coords_name[i], 'green')
        earth_plt = plt.Circle((0, 0, 0), R, color='g')
        ax.add_patch(earth_plt)
        ax.plot(*coords_name[i], 'red')
        ax.set(xlabel=f'{plots_2d[i][0]}-coordinate', ylabel=f'{plots_2d[i][1]}-coordinate', title=title + f'_{plots_2d[i][0]}{plots_2d[i][1]}')
        ax.grid()
        fig.tight_layout()

        if isinstance(name + f'_{plots_2d[i][0]}{plots_2d[i][1]}' + ext, str):
            fig.savefig(fname=name + f'_{plots_2d[i][0]}{plots_2d[i][1]}' + ext, bbox_inches='tight')



if __name__ == "__main__":
    plot_projectile(trajectory_path='satellite/trajectory/trajectory_1.csv', title=f"satellite", save_name='satellite/img/trajectory_1.png')