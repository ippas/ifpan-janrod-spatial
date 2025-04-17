import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re

def plot_boxplot_with_jitter(
    data,
    x,
    y,
    jitter_strength=0.1,
    point_size=6,
    sort_by=None,
    top_n_labels=5,
    show_top_genes=False,
    ax=None
):
    """
    Creates a boxplot with jittered data points colored by mouse genotype.
    Optionally adds gene labels for top N genes per category.

    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame containing the data to plot.
    x : str
        Column name to use for the X-axis (categorical).
    y : str
        Column name to use for the Y-axis (numeric).
    jitter_strength : float, optional
        The amount of jitter to apply to the data points (default is 0.1).
    point_size : float, optional
        Size of the jittered data points (default is 6).
    sort_by : str or None, optional
        Determines sorting of the X-axis categories. Options:
        "sample" (by sample number), "mean", "median", "max", "IQR" or None.
    top_n_labels : int, optional
        Number of top genes to label per X category (used if show_top_genes=True).
    show_top_genes : bool, optional
        Whether to display the top gene labels per X-axis category (default is False).
    ax : matplotlib.axes.Axes or None, optional
        Axis object to plot on. If None, a new figure and axis will be created.

    Returns:
    --------
    None
    """
    # Define custom color palette for genotypes
    genotype_palette = {"wtwt": "#1f77b4", "wtdel": "#d62728"}

    # Sort X-axis labels
    if sort_by == "sample":
        # Sort by sample number using regex
        order = sorted(data[x].unique(), key=lambda s: int(re.search(r'Nr(\d+)', s).group(1)))
        title_sort_text = " (sorted by sample Nr)"
    elif sort_by:
        # Calculate sort metric per category
        if sort_by == "mean":
            stat_values = data.groupby(x)[y].mean()
        elif sort_by == "median":
            stat_values = data.groupby(x)[y].median()
        elif sort_by == "max":
            stat_values = data.groupby(x)[y].max()
        elif sort_by == "IQR":
            stat_values = data.groupby(x)[y].apply(lambda x: x.quantile(0.75) - x.quantile(0.25))
        else:
            raise ValueError("Invalid value for sort_by.")
        # Sort categories by selected metric
        sorted_order = stat_values.sort_values(ascending=False).index
        data[x] = pd.Categorical(data[x], categories=sorted_order, ordered=True)
        order = sorted_order
        title_sort_text = f" (sorted by {sort_by})"
    else:
        # No sorting
        order = data[x].unique()
        title_sort_text = ""

    # Create new figure and axis if none provided
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 9))
        created_fig = True

    # Draw boxplot without outliers
    sns.boxplot(
        data=data,
        x=x,
        y=y,
        color="white",
        fliersize=0,
        linewidth=1.5,
        order=order,
        ax=ax
    )

    # Add jittered data points colored by genotype
    sns.stripplot(
        data=data,
        x=x,
        y=y,
        hue="mouse_genotype",
        palette=genotype_palette,
        dodge=False,
        jitter=jitter_strength,
        size=point_size,
        alpha=0.8,
        order=order,
        ax=ax
    )

    # Add gene labels above boxes for top N genes per category
    if show_top_genes:
        ymax = data[y].max()
        y_offset = ymax * 1.02
        ax.set_ylim(top=ymax * 1.25)  # Add space above plot
        for category in order:
            subset = data[data[x] == category].nlargest(top_n_labels, y)
            gene_list = subset['gene_symbol'].tolist()
            gene_text = '\n'.join(gene_list)
            x_pos = list(order).index(category)
            ax.text(
                x_pos - 0.4,
                y_offset,
                gene_text,
                ha='left',
                va='bottom',
                fontsize=6
            )

    # Axis labels and aesthetics
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel(y)
    ax.set_xlabel(x)
    ax.legend(title="Mouse genotype", bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_title(f"Max coverage per gene across {x} and {y}{title_sort_text}", fontsize=12, pad=15)

    # If this is a standalone plot, show the figure
    if created_fig:
        plt.tight_layout(rect=[0, 0, 1, 0.92])
        plt.suptitle(f"Max coverage per gene across {x} and {y}{title_sort_text}", fontsize=14, y=0.98)
        plt.show()

# # version 1
# def plot_boxplot_with_jitter(data, x, y, jitter_strength=0.1, point_size=6, sort_by=None, top_n_labels=5, show_top_genes=False):
#     genotype_palette = {"wtwt": "#1f77b4", "wtdel": "#d62728"}

#     if sort_by == "sample":
#         order = sorted(data[x].unique(), key=lambda s: int(re.search(r'Nr(\d+)', s).group(1)))
#         title_sort_text = " (sorted by sample Nr)"
#     elif sort_by:
#         if sort_by == "mean":
#             stat_values = data.groupby(x)[y].mean()
#         elif sort_by == "median":
#             stat_values = data.groupby(x)[y].median()
#         elif sort_by == "max":
#             stat_values = data.groupby(x)[y].max()
#         elif sort_by == "IQR":
#             stat_values = data.groupby(x)[y].apply(lambda x: x.quantile(0.75) - x.quantile(0.25))
#         else:
#             raise ValueError("Invalid value for sort_by.")

#         sorted_order = stat_values.sort_values(ascending=False).index
#         data[x] = pd.Categorical(data[x], categories=sorted_order, ordered=True)
#         order = sorted_order
#         title_sort_text = f" (sorted by {sort_by})"
#     else:
#         order = data[x].unique()
#         title_sort_text = ""

#     plt.figure(figsize=(16, 9))

#     sns.boxplot(
#         data=data, x=x, y=y, color="white",
#         fliersize=0, linewidth=1.5, order=order
#     )

#     ax = sns.stripplot(
#         data=data, x=x, y=y, hue="mouse_genotype",
#         palette=genotype_palette, dodge=False,
#         jitter=jitter_strength, size=point_size,
#         alpha=0.8, order=order
#     )

#     if show_top_genes:
#         ymax = data[y].max()
#         y_offset = ymax * 1.08

#         for category in order:
#             subset = data[data[x] == category].nlargest(top_n_labels, y)
#             gene_list = subset['gene_symbol'].tolist()
#             gene_text = '\n'.join(gene_list)
#             x_pos = list(order).index(category)
#             plt.text(
#                 x_pos - 0.4, y_offset, gene_text,
#                 ha='left', va='bottom', fontsize=7
#             )

#     plt.xticks(rotation=45, ha="right")
#     plt.ylabel(y)
#     plt.xlabel(x)

#     plt.legend(title="Mouse genotype", bbox_to_anchor=(1.05, 1), loc='upper left')
#     plt.tight_layout(rect=[0, 0, 1, 0.9])
#     plt.suptitle(f"Max coverage per gene across {x} and {y}{title_sort_text}", fontsize=14, y=0.92)
#     plt.show()

# # plot_boxplot_with_jitter(
# #     data,  # Your DataFrame
# #     x="label",  # Column for X-axis
# #     y="norm_by_gapdh",  # Column for Y-axis
# #     jitter_strength=0.2,  # Adjust jitter spread (larger values = more spread)
# #     point_size=3,  # Adjust the size of the jittered points
# #     sort_by="sample",  # Sort the boxes by the IQR of max_coverage in descending order
# #     top_n_labels=5,  # Show top 5 genes for each sample
# #     show_top_genes=True  # Show top 5 genes for each sample
# # )

# # version 2
# def plot_boxplot_with_jitter(
#     data,
#     x,
#     y,
#     jitter_strength=0.1,
#     point_size=6,
#     sort_by=None,
#     top_n_labels=5,
#     show_top_genes=False,
#     ax=None
# ):
#     genotype_palette = {"wtwt": "#1f77b4", "wtdel": "#d62728"}

#     # Sortowanie etykiet na osi X
#     if sort_by == "sample":
#         order = sorted(data[x].unique(), key=lambda s: int(re.search(r'Nr(\d+)', s).group(1)))
#         title_sort_text = " (sorted by sample Nr)"
#     elif sort_by:
#         if sort_by == "mean":
#             stat_values = data.groupby(x)[y].mean()
#         elif sort_by == "median":
#             stat_values = data.groupby(x)[y].median()
#         elif sort_by == "max":
#             stat_values = data.groupby(x)[y].max()
#         elif sort_by == "IQR":
#             stat_values = data.groupby(x)[y].apply(lambda x: x.quantile(0.75) - x.quantile(0.25))
#         else:
#             raise ValueError("Invalid value for sort_by.")
#         sorted_order = stat_values.sort_values(ascending=False).index
#         data[x] = pd.Categorical(data[x], categories=sorted_order, ordered=True)
#         order = sorted_order
#         title_sort_text = f" (sorted by {sort_by})"
#     else:
#         order = data[x].unique()
#         title_sort_text = ""

#     # Jeśli nie podano ax — twórz nową figurę i oś
#     created_fig = False
#     if ax is None:
#         fig, ax = plt.subplots(figsize=(16, 9))
#         created_fig = True

#     # Boxplot
#     sns.boxplot(
#         data=data,
#         x=x,
#         y=y,
#         color="white",
#         fliersize=0,
#         linewidth=1.5,
#         order=order,
#         ax=ax
#     )

#     # Jittered punkty
#     sns.stripplot(
#         data=data,
#         x=x,
#         y=y,
#         hue="mouse_genotype",
#         palette=genotype_palette,
#         dodge=False,
#         jitter=jitter_strength,
#         size=point_size,
#         alpha=0.8,
#         order=order,
#         ax=ax
#     )

#     # Górne etykiety z nazwami genów
#     if show_top_genes:
#         ymax = data[y].max()
#         y_offset = ymax * 1.02
#         ax.set_ylim(top=ymax * 1.25)  # Dodaj przestrzeń nad wykresem
#         for category in order:
#             subset = data[data[x] == category].nlargest(top_n_labels, y)
#             gene_list = subset['gene_symbol'].tolist()
#             gene_text = '\n'.join(gene_list)
#             x_pos = list(order).index(category)
#             ax.text(
#                 x_pos - 0.4,
#                 y_offset,
#                 gene_text,
#                 ha='left',
#                 va='bottom',
#                 fontsize=6
#             )

#     # Estetyka
#     ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
#     ax.set_ylabel(y)
#     ax.set_xlabel(x)
#     ax.legend(title="Mouse genotype", bbox_to_anchor=(1.05, 1), loc='upper left')
#     ax.set_title(f"Max coverage per gene across {x} and {y}{title_sort_text}", fontsize=12, pad=15)

#     # Jeśli to wykres samodzielny — dodaj suptitle i pokaż
#     if created_fig:
#         plt.tight_layout(rect=[0, 0, 1, 0.92])
#         plt.suptitle(f"Max coverage per gene across {x} and {y}{title_sort_text}", fontsize=14, y=0.98)
#         plt.show()
