#!/usr/bin/env python3
"""
Simple script for managing command-line input.
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description="Manage command-line input")

    parser.add_argument("input_file", help="Input file path")
    parser.add_argument("-o", "--output", default="output.txt", help="Output file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")

    args = parser.parse_args()

    print(f"Input: {args.input_file}")
    print(f"Output: {args.output}")
    if args.verbose:
        print("Verbose mode enabled")


if __name__ == "__main__":
    main()
