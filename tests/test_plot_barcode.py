import pytest
import pandas as pd
from pathlib import Path
from barcodeforge.plot_barcode import plot_barcode_altair, create_barcode_plot
import altair as alt
import json


@pytest.fixture
def sample_barcode_df_long():
    data = {
        "Lineage": ["L1", "L1", "L2", "L2", "L3", "L3"],
        "Mutation": ["M1", "M2", "M1", "M3", "M2", "M3"],
        "z": [1, 0, 1, 1, 0, 1],
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_barcode_csv_file(tmp_path):
    file_path = tmp_path / "sample_barcodes.csv"
    data = {
        "M1": {"L1": 1, "L2": 1, "L3": 0},
        "M2": {"L1": 0, "L2": 0, "L3": 1},  # Changed L3 M2 to 1 to have a present value
        "M3": {"L1": 0, "L2": 1, "L3": 1},
    }
    df = pd.DataFrame(data)
    df.index.name = "Lineage"
    df.to_csv(file_path)
    return file_path


def test_plot_barcode_altair_creates_html(sample_barcode_df_long, tmp_path):
    output_html = tmp_path / "test_plot.html"
    plot_barcode_altair(sample_barcode_df_long, str(output_html))
    assert output_html.exists()
    # Further checks could involve parsing the HTML or checking its content,
    # but that can be complex. For now, just check creation.
    with open(output_html, "r") as f:
        content = f.read()
        assert '<div id="vis">' in content  # Basic check for Altair plot presence


def test_plot_barcode_altair_chart_structure(sample_barcode_df_long, tmp_path, mocker):
    # Mock the save method on the VConcatChart object that vconcat returns
    mock_vconcat_save = mocker.patch("altair.VConcatChart.save")

    # Spy on alt.Chart to check its calls without altering its behavior for chart creation
    spy_chart = mocker.spy(alt, "Chart")

    output_html = tmp_path / "test_plot_structure.html"
    plot_barcode_altair(sample_barcode_df_long, str(output_html))

    # Check if the 'save' method on VConcatChart was called
    mock_vconcat_save.assert_called_once_with(str(output_html))

    # Check that alt.Chart was called twice (for box1 and box2)
    assert spy_chart.call_count == 2

    # Check the first call to alt.Chart (for box1, then box2 due to vconcat order)
    # The actual order in vconcat is box2, box1. So spy_chart.call_args_list[0] is box2, spy_chart.call_args_list[1] is box1
    pd.testing.assert_frame_equal(
        spy_chart.call_args_list[0][0][0], sample_barcode_df_long
    )  # args[0] is the df
    pd.testing.assert_frame_equal(
        spy_chart.call_args_list[1][0][0], sample_barcode_df_long
    )  # args[0] is the df


def test_create_barcode_plot_integration(sample_barcode_csv_file, tmp_path, mocker):
    output_html = tmp_path / "integrated_plot.html"

    # Mock the underlying plot_barcode_altair to simplify this integration test
    # and focus on the file reading and data transformation part of create_barcode_plot
    mock_plot_altair = mocker.patch("barcodeforge.plot_barcode.plot_barcode_altair")

    create_barcode_plot(str(sample_barcode_csv_file), str(output_html))

    # Check that plot_barcode_altair was called with the correct arguments
    mock_plot_altair.assert_called_once()
    call_args = mock_plot_altair.call_args[0]
    df_long_arg = call_args[0]
    output_path_arg = call_args[1]

    assert isinstance(df_long_arg, pd.DataFrame)
    assert "Lineage" in df_long_arg.columns
    assert "Mutation" in df_long_arg.columns
    assert "z" in df_long_arg.columns
    assert output_path_arg == str(output_html)

    # Verify the transformation from wide to long format
    # L1, M1 should be 1
    assert (
        df_long_arg[
            (df_long_arg["Lineage"] == "L1") & (df_long_arg["Mutation"] == "M1")
        ]["z"].iloc[0]
        == 1
    )
    # L1, M2 should be 0
    assert (
        df_long_arg[
            (df_long_arg["Lineage"] == "L1") & (df_long_arg["Mutation"] == "M2")
        ]["z"].iloc[0]
        == 0
    )


def test_create_barcode_plot_file_creation(sample_barcode_csv_file, tmp_path):
    # This is a more direct test of file creation without mocking the altair plot function itself
    output_html = tmp_path / "direct_integrated_plot.html"
    create_barcode_plot(str(sample_barcode_csv_file), str(output_html))
    assert output_html.exists()
    with open(output_html, "r") as f:
        content = f.read()
        assert '<div id="vis">' in content
        # Check for some data from the plot to ensure it's not just an empty template
        assert "L1" in content
        assert "M1" in content
