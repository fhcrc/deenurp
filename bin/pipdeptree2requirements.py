#!/usr/bin/env python

"""Generate a requirements.txt file containing frozen dependencies
with secondary (or higher-order) dependencies first.

Usage:

  pip install pipdeptree
  pipdeptree -f | bin/pipdeptree2requirements.py

"""

import sys
import re

lines = []
for line in sys.stdin:
    spaces = re.findall(r'^[ ]+', line)
    lines.append((-1 * len(spaces[0]) if spaces else 0, line.strip()))

deps = set()
for __, pkg in sorted(lines):
    if pkg not in deps:
        print pkg
        deps.add(pkg)
