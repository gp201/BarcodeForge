import pytest
import click
import subprocess
from unittest.mock import MagicMock, patch, call
from rich.console import Console
from barcodeforge.utils import (
    sortFun,
    resolve_tree_format,
    run_subprocess_command,
    STYLES,
)


def test_sortFun():
    assert sortFun("A123B") == 123
    assert sortFun("C456D") == 456
    assert sortFun("G789E") == 789


def test_resolve_tree_format_specified_newick():
    assert resolve_tree_format("any.tree", "newick", None, False) == "newick"


def test_resolve_tree_format_specified_nexus():
    assert resolve_tree_format("any.tree", "nexus", None, False) == "nexus"


def test_resolve_tree_format_infer_nwk():
    assert resolve_tree_format("test.nwk", None, None, False) == "newick"


def test_resolve_tree_format_infer_newick():
    assert resolve_tree_format("test.newick", None, None, False) == "newick"


def test_resolve_tree_format_infer_nexus():
    assert resolve_tree_format("test.nexus", None, None, False) == "nexus"


def test_resolve_tree_format_unknown_extension():
    with pytest.raises(click.Abort):
        resolve_tree_format("test.txt", None, None, False)


def test_resolve_tree_format_debug_output():
    mock_console = MagicMock(spec=Console)
    resolve_tree_format("some.nwk", None, mock_console, debug=True)
    mock_console.print.assert_any_call(
        f"[{STYLES['warning']}]Resolved tree format for 'some.nwk': newick[/]"
    )


@patch("subprocess.run")
def test_run_subprocess_command_success(mock_subproc_run):
    mock_subproc_run.return_value = subprocess.CompletedProcess(
        args=["test_cmd"], returncode=0, stdout="Success output", stderr=""
    )
    mock_console = MagicMock(spec=Console)
    result = run_subprocess_command(
        ["test_cmd"],
        mock_console,
        debug=False,
        success_message="Command executed successfully",
    )
    assert result is True
    mock_console.print.assert_called_once_with(
        f"[{STYLES['success']}]Command executed successfully[/]"
    )


@patch("subprocess.run")
def test_run_subprocess_command_failure_called_process_error(mock_subproc_run):
    # This test covers cases where subprocess.run(check=True) would raise CalledProcessError
    mock_subproc_run.side_effect = subprocess.CalledProcessError(
        returncode=1, cmd=["fail_cmd_cpe"], output="out", stderr="Error output cpe"
    )
    mock_console = MagicMock(spec=Console)
    with pytest.raises(click.Abort):
        run_subprocess_command(
            ["fail_cmd_cpe"],
            mock_console,
            debug=False,
            error_message_prefix="Test error CPE",
        )

    expected_calls = [
        call(
            f"[{STYLES['error']}]Test error CPE fail_cmd_cpe: Command '['fail_cmd_cpe']' returned non-zero exit status 1.[/]"
        ),
        call(
            f"[{STYLES['dim']}]fail_cmd_cpe stdout:\\\\nout[/]"
        ),  # Corrected: output is stdout
        call(f"[{STYLES['dim']}]fail_cmd_cpe stderr:\\\\nError output cpe[/]"),
    ]
    mock_console.print.assert_has_calls(expected_calls)


@patch("subprocess.run")
def test_run_subprocess_command_file_not_found(mock_subproc_run):
    mock_subproc_run.side_effect = FileNotFoundError("Command not found")
    mock_console = MagicMock(spec=Console)
    with pytest.raises(click.Abort):
        run_subprocess_command(
            ["non_existent_cmd"],
            mock_console,
            debug=False,
            error_message_prefix="FNF error",
        )

    mock_console.print.assert_called_once_with(
        f"[{STYLES['error']}]FNF error: non_existent_cmd command not found. Please ensure it is installed and in your PATH.[/]"
    )


@patch("subprocess.run")
def test_run_subprocess_command_success_debug(mock_subproc_run):
    mock_subproc_run.return_value = subprocess.CompletedProcess(
        args=["debug_cmd_success"],
        returncode=0,
        stdout="Debug success output",
        stderr="Debug success stderr",
    )
    mock_console = MagicMock(spec=Console)
    result = run_subprocess_command(
        ["debug_cmd_success", "arg1"],
        mock_console,
        debug=True,
        success_message="Debug success",
    )
    assert result is True
    expected_calls = [
        call(f"[{STYLES['debug']}]Running command: debug_cmd_success arg1[/]"),
        call(f"[{STYLES['dim']}]debug_cmd_success stdout:\\\\nDebug success output[/]"),
        call(f"[{STYLES['dim']}]debug_cmd_success stderr:\\\\nDebug success stderr[/]"),
        call(f"[{STYLES['success']}]Debug success[/]"),
    ]
    mock_console.print.assert_has_calls(expected_calls, any_order=False)


@patch("subprocess.run")
def test_run_subprocess_command_success_debug_empty_stderr(mock_subproc_run):
    mock_subproc_run.return_value = subprocess.CompletedProcess(
        args=["debug_cmd_success_no_stderr"],
        returncode=0,
        stdout="Debug success output",
        stderr="",
    )
    mock_console = MagicMock(spec=Console)
    result = run_subprocess_command(
        ["debug_cmd_success_no_stderr"],
        mock_console,
        debug=True,
        success_message="Debug success no stderr",
    )
    assert result is True
    expected_calls = [
        call(f"[{STYLES['debug']}]Running command: debug_cmd_success_no_stderr[/]"),
        call(
            f"[{STYLES['dim']}]debug_cmd_success_no_stderr stdout:\\\\nDebug success output[/]"
        ),
        # No call for stderr as it's empty
        call(f"[{STYLES['success']}]Debug success no stderr[/]"),
    ]
    mock_console.print.assert_has_calls(expected_calls, any_order=False)
    # Verify stderr was not printed.
    # We check that if a call contains "stderr:\\n", it must be followed by another character,
    # meaning it's not an empty stderr print.
    for acall in mock_console.print.call_args_list:
        call_str = acall[0][0]
        # This assertion is a bit tricky. We want to ensure that if "stderr:\\n" is present,
        # it's not *just* "stderr:\\n" (or "stderr:\\n[/]" with a style closing tag).
        # It should have content after "stderr:\\n".
        # A simpler way is to count calls that would match an empty stderr print.
        # However, the current structure of the code in utils.py ensures that
        # if stderr is empty, the print call for stderr is skipped entirely.
        # So, we just need to ensure no call looks like an empty stderr print.
        assert not (
            f" stderr:\\\\n</]" in call_str and call_str.endswith(f" stderr:\\\\n</]")
        )
        assert not (
            f" stderr:\\\\n[/{STYLES['dim']}]"
            == call_str[-len(f" stderr:\\\\n[/{STYLES['dim']}]") :]
        )  # check end of string


@patch("subprocess.run")
def test_run_subprocess_command_failure_debug(mock_subproc_run):
    cmd_list = ["debug_cmd_fail", "arg_fail"]
    mock_subproc_run.side_effect = subprocess.CalledProcessError(
        returncode=1,
        cmd=cmd_list,
        output="Debug fail stdout",
        stderr="Debug fail stderr",
    )
    mock_console = MagicMock(spec=Console)
    with pytest.raises(click.Abort):
        run_subprocess_command(
            cmd_list, mock_console, debug=True, error_message_prefix="Debug fail error"
        )

    expected_calls_in_order = [
        call(f"[{STYLES['debug']}]Running command: {' '.join(cmd_list)}[/]"),
        call(
            f"[{STYLES['error']}]Debug fail error {cmd_list[0]}: Command '{cmd_list}' returned non-zero exit status 1.[/]"
        ),
        call(f"[{STYLES['dim']}]{cmd_list[0]} stdout:\\\\nDebug fail stdout[/]"),
        call(f"[{STYLES['dim']}]{cmd_list[0]} stderr:\\\\nDebug fail stderr[/]"),
    ]
    mock_console.print.assert_has_calls(expected_calls_in_order, any_order=False)


@patch("subprocess.run")
def test_run_subprocess_command_no_console_success(mock_subproc_run):
    mock_subproc_run.return_value = subprocess.CompletedProcess(
        args=["test_cmd_no_console"], returncode=0, stdout="Success output", stderr=""
    )
    # Pass None for console
    result = run_subprocess_command(
        ["test_cmd_no_console"],
        None,
        debug=False,
        success_message="Command executed successfully",
    )
    assert result is True
    # mock_console.print should not have been called, as console is None.
    # No direct assertion for no calls on None, but the command should complete without AttributeError.


@patch("subprocess.run")
def test_run_subprocess_command_no_console_failure(mock_subproc_run):
    mock_subproc_run.side_effect = subprocess.CalledProcessError(
        returncode=1, cmd=["fail_cmd_no_console"], output="", stderr="Error output"
    )
    with pytest.raises(click.Abort):
        # Pass None for console
        run_subprocess_command(
            ["fail_cmd_no_console"],
            None,
            debug=False,
            error_message_prefix="Test error no console",
        )
    # No direct assertion for no calls on None, but the command should raise click.Abort without AttributeError.


@patch("subprocess.run")
def test_run_subprocess_command_no_console_debug_success(mock_subproc_run):
    mock_subproc_run.return_value = subprocess.CompletedProcess(
        args=["debug_cmd_success_no_console"],
        returncode=0,
        stdout="Debug success output",
        stderr="",
    )
    # Pass None for console, debug True
    result = run_subprocess_command(
        ["debug_cmd_success_no_console"],
        None,
        debug=True,
        success_message="Debug success",
    )
    assert result is True
    # Debug messages and success messages are conditional on console being non-None.
    # Primary check is that it runs without error.
