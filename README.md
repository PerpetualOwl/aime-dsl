# GeometryDSL Project

## Overview

GeometryDSL is a domain-specific language (DSL) designed to facilitate the creation, manipulation, and solving of geometric problems. It leverages symbolic mathematics to define geometric entities, apply constraints, and compute solutions using Python's `sympy` library.

## Features

- **Define Geometric Entities**: Easily specify points, lines, circles, and angles.
- **Constraints Handling**: Set and solve constraints on distances, angles, and other geometric properties.
- **Automatic Solving**: Calculate unknowns using symbolic algebra and geometry principles.
- **Flexible Use Cases**: Suitable for educational demonstrations, algorithm testing, and geometric problem-solving.

## Installation

To use the GeometryDSL, ensure you have Python and the necessary dependencies installed:

```sh
pip install sympy
```

## Usage

Start by defining geometric entities and constraints using the `GeometryDSL` class. Below is a basic example demonstrating how to define a triangle and calculate its properties:

```python
from sympy import *
from sympy.geometry import *

# Instantiate the DSL
dsl = GeometryDSL()

# Define points
dsl.define_point("A")
dsl.define_point("B")
dsl.define_point("C")

# Define lines
dsl.define_line("AB", "A", "B")
dsl.define_line("BC", "B", "C")
dsl.define_line("AC", "A", "C")

# Set lengths
dsl.set_length("AB", 5)
dsl.set_length("BC", 7)
dsl.set_length("AC", 8)

# Solve for angles
angles = dsl.solve("BAC")
print("Angle BAC:", angles)
```

## LLM Usage

A system prompt is provided in `system_prompt.txt`. It is suggested to use this in conjunction with providing the python declarations of dataclasses in a query to an LLM to generate code that is usable with the DSL. Examples 2 and 3 below were generated with the help of an LLM.

## Examples

- **example1**: Solving a basic triangle problem, calculating coordinates and angles.
- **example2**: Working with circles and determining if a point lies inside it.
- **example3**: Define an equilateral triangle and solve for coordinates.