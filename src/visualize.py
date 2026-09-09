import matplotlib.pyplot as plt
import seaborn as sns
import os 


def create_visualizations(df):
    """
    Generates, formats, and saves a summary visualization dashboard 
    containing GC percentage bar plots and length vs. molecular weight scatter plots.

    Args:
        df (pd.DataFrame): The DataFrame containing processed sequence statistics.

    Returns:
        None
    """

    sns.set_theme(style = "whitegrid")

    fig, axes = plt.subplots(1, 2, figsize = (14, 6))

    sns.barplot( 
        data = df, 
        x= "id",
        y = "GC_Percent", 
        ax= axes[0],
        palette= "viridis",
        legend = False
    )

    axes[0].set_title("GC Percentage per Sequence")
    axes[0].set_ylabel("GC Content (%)")

    axes[0].set_xticks(range(len(df["id"])))
    axes[0].set_xticklabels(df["id"], rotation=45, ha="right")

    sns.scatterplot(
        data=df, 
        x="length", 
        y="Mol_Weight_(Da)", 
        ax=axes[1], 
        hue=df["GC_Percent"], 
        palette="magma", 
        s=100
    )

    axes[1].set_title("Sequence Length vs. Molecular Weight")
    axes[1].set_xlabel("Sequence Length (bp)")
    axes[1].set_ylabel("Molecular Weight (Da)")
    
    plt.tight_layout()
    os.makedirs("outputs", exist_ok=True)
    plot_filename = "outputs/sequence_analysis_dashboard.png"
    plt.savefig(plot_filename)
    print(f"\nDashboard saved as '{plot_filename}'")
    plt.close()


