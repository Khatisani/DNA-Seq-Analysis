import matplotlib.pyplot as plt
import seaborn as sns

#generates and saves a matplotlib/seaborn visualization from the DataFrame
def create_visualizations(df):
    sns.set_theme(style = "whitegrid")

    fig, axes = plt.subplots(1, 2, fiqsize = (14, 6))

    sns.barplot( 
        data = df, 
        x= "id",
        y = "gc%", 
        ax= axes[0],
        palette= "viridis"
    )

    axes[0].set_title("GC Percentage per Sequence")
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha="right")
    axes[0].set_ylabel("GC Content (%)")

    sns.scatterplot(
        data=df, 
        x="length", 
        y="Mol_Weight_(Da)", 
        ax=axes[1], 
        hue="GC%", 
        palette="magma", 
        s=100
    )

    axes[1].set_title("Sequence Length vs. Molecular Weight")
    axes[1].set_xlabel("Sequence Length (bp)")
    axes[1].set_ylabel("Molecular Weight (Da)")
    
    plt.tight_layout()
    plot_filename = "sequence_analysis_dashboard.png"
    plt.savefig(plot_filename)
    print(f"\nDashboard saved as '{plot_filename}'")
    plt.close()
    

