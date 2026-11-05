#!/usr/bin/env python3
"""
General utility template for constructing input arguments in Python.

This module can be:
  1. Imported by other scripts to reuse argument handling patterns
  2. Run as a standalone script for testing

Usage as module:
  from input_argument_template import add_file_arguments, add_format_arguments

Usage as script:
  python input_argument_template.py -i input.txt -o output.txt --format json
"""

import argparse
from pathlib import Path
from typing import Optional


def add_file_arguments(parser: argparse.ArgumentParser,
                       input_required: bool = True) -> None:
    """Add common file-related arguments to parser.

    Args:
        parser: ArgumentParser instance to add arguments to
        input_required: Whether input file is required (positional vs optional)
    """
    if input_required:
        parser.add_argument(
            "input_file",
            type=Path,
            help="Path to input file"
        )
    else:
        parser.add_argument(
            "-i", "--input",
            type=Path,
            required=False,
            help="Path to input file"
        )

    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("output.txt"),
        help="Path to output file (default: output.txt)"
    )


def add_format_arguments(parser: argparse.ArgumentParser,
                        formats: Optional[list] = None) -> None:
    """Add common format-related arguments to parser.

    Args:
        parser: ArgumentParser instance to add arguments to
        formats: List of valid format choices (default: json, csv, txt)
    """
    if formats is None:
        formats = ["json", "csv", "txt"]

    parser.add_argument(
        "--format",
        choices=formats,
        default=formats[0],
        help=f"Output format (default: {formats[0]})"
    )


def add_verbosity_arguments(parser: argparse.ArgumentParser) -> None:
    """Add common verbosity arguments to parser.

    Args:
        parser: ArgumentParser instance to add arguments to
    """
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )

    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress output"
    )


def create_parser(description: str = "General utility script",
                  add_defaults: bool = True) -> argparse.ArgumentParser:
    """Create and configure argument parser.

    Args:
        description: Description of the script
        add_defaults: Whether to add default common arguments

    Returns:
        Configured ArgumentParser instance
    """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    if add_defaults:
        add_file_arguments(parser)
        add_format_arguments(parser)
        add_verbosity_arguments(parser)

    return parser


def process_file(input_file: Path, output_file: Path,
                output_format: str, verbose: bool = False) -> None:
    """Process the input file.

    Args:
        input_file: Path to input file
        output_file: Path to output file
        output_format: Output format
        verbose: Whether to print verbose output
    """
    if verbose:
        print(f"Input: {input_file}")
        print(f"Output: {output_file}")
        print(f"Format: {output_format}")

    # Implementation goes here
    pass


def main() -> None:
    """Main entry point when run as standalone script."""
    parser = create_parser(description="Template utility script")
    args = parser.parse_args()

    input_file = getattr(args, 'input_file', None) or args.input
    output_file = args.output
    output_format = args.format
    verbose = args.verbose and not args.quiet

    process_file(input_file, output_file, output_format, verbose)


if __name__ == "__main__":
    main()
