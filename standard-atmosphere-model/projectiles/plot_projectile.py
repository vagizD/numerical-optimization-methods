import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def plot_projectile(trajectory_path: str, title: str, save_name: str = None):
    points = pd.read_csv(trajectory_path)

    x0, y0 = points.x.values[0], points.y.values[0]
    x_end, y_end = points.x.values[-1], points.y.values[-1]
    y_min_idx, y_max_idx = points.y.values.argmin(), points.y.values.argmax()

    fig, ax = plt.subplots(figsize=(10, 6))

    sns.scatterplot(x=points.x, y=points.y, label="Projectile", ax=ax)
    ax.grid()

    # ground
    ax.axhline(y=y0, xmin=x0, xmax=x_end, c='k')
    # range
    ax.axvline(x=points.x.values[y_max_idx], ymin=y0, ymax=points.y.values[y_max_idx],
               label=f"Range = {points.y.values[y_max_idx]:.3f}", ls='--', c='r')

    ax.set(xlabel='x-coordinate', ylabel='y-coordinate', title=title)
    fig.tight_layout()

    if isinstance(save_name, str):
        fig.savefig(fname=save_name, dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    from os.path import join

    save_name = join('img', 'trajectory1.png')
    trajectory_path = join('trajectory', 'trajectory1.csv')

    plot_projectile(trajectory_path=trajectory_path, title="First Attempt", save_name=save_name)
