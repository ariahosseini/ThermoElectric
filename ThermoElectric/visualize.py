import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Set up Seaborn visualization style
sns.set()
sns.set_style("white", {"xtick.major.size": 2, "ytick.major.size": 2})
sns.set_context("paper", font_scale=2, rc={"lines.linewidth": 4})

def visualize(array_to_plot: np.ndarray, linestyle: str, marker: str, pic_name: str,
              x_label: str = None, y_label: str = None) -> None:
    """
    Plot a two-dimensional numpy array with specified style and save it as an image.

    Parameters
    ----------
    array_to_plot : np.ndarray
        A 2D numpy array where the first row is x-values and the second row is y-values.
    linestyle : str
        Line style for the plot (e.g., '-', '--').
    marker : str
        Marker style for the plot (e.g., 'o', '^').
    pic_name : str
        Name of the file where the plot will be saved (should include file extension, e.g., 'plot.png').
    x_label : str, optional
        Label for the x-axis (default is None).
    y_label : str, optional
        Label for the y-axis (default is None).

    Returns
    -------
    None
    """

    # Create directory if it doesn't exist
    dir_name = 'Figs'
    os.makedirs(dir_name, exist_ok=True)

    # Plot settings
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.grid(False)  # Disable the grid for a cleaner plot
    ax.plot(array_to_plot[0], array_to_plot[1], linestyle=linestyle, marker=marker, color='maroon',
            markersize=6, linewidth=2, markerfacecolor='white', markeredgecolor='maroon', markeredgewidth=1)

    # Axis labels and formatting
    if x_label:
        ax.set_xlabel(x_label, fontsize=20)
    if y_label:
        ax.set_ylabel(y_label, fontsize=20, labelpad=15)
    
    ax.tick_params(axis="x", labelsize=20)
    ax.tick_params(axis="y", labelsize=20)
    
    # Scientific notation for large/small numbers
    plt.ticklabel_format(axis="both", style="sci")

    # Save figure
    fig.tight_layout()
    plt.savefig(os.path.join(dir_name, pic_name), dpi=100)
    plt.close()
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Set up Seaborn visualization style
sns.set()
sns.set_style("white", {"xtick.major.size": 2, "ytick.major.size": 2})
sns.set_context("paper", font_scale=2, rc={"lines.linewidth": 4})

def visualize(array_to_plot: np.ndarray, linestyle: str, marker: str, pic_name: str,
              x_label: str = None, y_label: str = None) -> None:
    """
    Plot a two-dimensional numpy array with specified style and save it as an image.

    Parameters
    ----------
    array_to_plot : np.ndarray
        A 2D numpy array where the first row is x-values and the second row is y-values.
    linestyle : str
        Line style for the plot (e.g., '-', '--').
    marker : str
        Marker style for the plot (e.g., 'o', '^').
    pic_name : str
        Name of the file where the plot will be saved (should include file extension, e.g., 'plot.png').
    x_label : str, optional
        Label for the x-axis (default is None).
    y_label : str, optional
        Label for the y-axis (default is None).

    Returns
    -------
    None
    """

    # Create directory if it doesn't exist
    dir_name = 'Figs'
    os.makedirs(dir_name, exist_ok=True)

    # Plot settings
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.grid(False)  # Disable the grid for a cleaner plot
    ax.plot(array_to_plot[0], array_to_plot[1], linestyle=linestyle, marker=marker, color='maroon',
            markersize=6, linewidth=2, markerfacecolor='white', markeredgecolor='maroon', markeredgewidth=1)

    # Axis labels and formatting
    if x_label:
        ax.set_xlabel(x_label, fontsize=20)
    if y_label:
        ax.set_ylabel(y_label, fontsize=20, labelpad=15)
    
    ax.tick_params(axis="x", labelsize=20)
    ax.tick_params(axis="y", labelsize=20)
    
    # Scientific notation for large/small numbers
    plt.ticklabel_format(axis="both", style="sci")

    # Save figure
    fig.tight_layout()
    plt.savefig(os.path.join(dir_name, pic_name), dpi=100)
    plt.close()
