"""Plot barcode from CSV file."""

import altair as alt
import pandas as pd
import math
from rich.console import Console
from .utils import sortFun, STYLES

console = Console()


def plot_barcode_altair(barcode_df_long: pd.DataFrame, output_html_path: str) -> None:
    """
    Generates and saves an Altair plot from a long-format barcode DataFrame.

    Args:
        barcode_df_long: Pandas DataFrame in long format with columns ["Lineage", "Mutation", "z"].
        output_html_path: Path to save the generated HTML plot.
    """
    # Map 0/1 to labels for discrete presence/absence
    barcode_df_long["Presence"] = barcode_df_long["z"].map({0: "Absent", 1: "Present"})
    # Determine dynamic dimensions based on unique categories in the dataframe
    n_lineages = len(barcode_df_long["Lineage"].unique())
    n_mutations = len(barcode_df_long["Mutation"].unique())

    # Calculate width and height dynamically
    # Adjusted calculations to prevent zero or excessively small dimensions
    width1 = max(
        200, round(math.log1p(n_mutations)) * 100
    )  # log1p for stability if n_mutations is small
    height1 = max(150, round(math.log1p(n_lineages)) * 60)

    width2 = max(300, n_mutations * 20)  # Ensure a minimum width
    height2 = max(200, n_lineages * 20)  # Ensure a minimum height

    console.print(
        f"[{STYLES['info']}]Plot dimensions: Zoomed-out ({width1}x{height1}), Zoomed-in ({width2}x{height2})[/{STYLES['info']}]"
    )

    box1 = (
        alt.Chart(barcode_df_long)
        .mark_rect(
            stroke="#FFFFFF",
            strokeWidth=0,
        )
        .encode(
            y="Lineage:O",
            x=alt.X("Mutation:O", sort=None, axis=alt.Axis(labels=False, tickSize=0)),
            color=alt.Color(
                "Presence:N",
                scale=alt.Scale(
                    domain=["Absent", "Present"], range=["#FFFFFF", "#000000"]
                ),
                legend=alt.Legend(title="Mutation status"),
            ),
            tooltip=["Lineage", "Mutation", "Presence"],
        )
        .properties(
            title=alt.TitleParams(
                text="Zoomed out Barcode", anchor="start", fontSize=20
            ),
            width=width1,
            height=height1,
        )
        .interactive()
    )

    box2 = (
        alt.Chart(barcode_df_long)
        .mark_rect(stroke="#BBB", strokeWidth=0.25)
        .encode(
            y="Lineage:O",
            x=alt.X("Mutation:O", sort=None, axis=alt.Axis(labels=True, labelAngle=45)),
            color=alt.Color(
                "Presence:N",
                scale=alt.Scale(
                    domain=["Absent", "Present"], range=["#FFFFFF", "#000000"]
                ),
                legend=alt.Legend(title="Mutation status"),
            ),
            tooltip=["Lineage", "Mutation", "Presence"],
        )
        .properties(
            title=alt.TitleParams(
                text="Zoomed in Barcode", anchor="start", fontSize=20
            ),
            width=width2,
            height=height2,
        )
        .interactive()
    )

    # Combine the two charts vertically
    combined_chart = alt.vconcat(box2, box1, spacing=10).resolve_scale(
        color="independent"
    )

    # Save the chart to an HTML file
    combined_chart.save(output_html_path)
    console.print(
        f"[{STYLES['success']}]Barcode plot saved to: {output_html_path}[/{STYLES['success']}]"
    )


def create_barcode_plot(input_file_path: str, output_file_path: str) -> None:
    """
    Reads a barcode CSV file, transforms it to long format, and generates an HTML plot.

    Args:
        input_file_path: Path to the input barcode CSV file.
        output_file_path: Path to save the generated HTML plot.
    """
    console.print(
        f"[{STYLES['info']}]Reading barcode data from: {input_file_path}[/{STYLES['info']}]"
    )
    barcode_df = pd.read_csv(input_file_path, header=0, index_col=0)

    # convert barcode dataframe to long format
    console.print(
        f"[{STYLES['info']}]Transforming barcode data to long format...[/{STYLES['info']}]"
    )
    barcode_df_long = barcode_df.stack().reset_index()
    barcode_df_long.columns = ["Lineage", "Mutation", "z"]

    plot_barcode_altair(barcode_df_long, output_file_path)
