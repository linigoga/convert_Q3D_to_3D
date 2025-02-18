# Quasi-3D to Full 3D Data Conversion Tool

This tool converts data files from quasi-3D simulations (like those produced by Osiris) to full 3D Cartesian data. It supports converting charge density, field, and raw particle data.  The conversion is performed using Cython for performance.

## Table of Contents

1.  [Introduction](#introduction)
2.  [Dependencies](#dependencies)
3.  [Installation](#installation)
4.  [Usage](#usage)
    *   [Command-Line Arguments](#command-line-arguments)
    *   [Examples](#examples)
5.  [Code Structure](#code-structure)
    *   `transform_data.py`
    *   `utils.pyx`
    *   `setup.py`
6.  [Workflow](#workflow)
7.  [Important Notes](#important-notes)
8.  [Troubleshooting](#troubleshooting)
9. [To do](#to-do)

## 1. Introduction

This tool takes simulation data in a cylindrical coordinate system (r, z, θ) with azimuthal symmetry (often used in quasi-3D simulations) and interpolates it to a full 3D Cartesian grid (x, y, z). This is useful for visualizing and analyzing data in standard 3D visualization tools. The tool supports multiple data types and provides options for specifying the output domain, timestep, and other parameters. The conversion is heavily optimized using Cython, providing significant speed improvements over a pure Python implementation.

## 2. Dependencies

You *must* have the following Python packages installed:

*   **numpy:** For numerical operations and array handling.
*   **scipy:** For interpolation (specifically `scipy.interpolate.interp1d` and `scipy.ndimage.map_coordinates`).
*   **h5py:** For reading and writing HDF5 files.
*   **Cython:** For compiling and using the optimized Cython code.
*   **setuptools:** For building the Cython extension.

It is **highly recommended** to use a virtual environment to manage these dependencies.

## 3. Installation

1.  **Clone the Repository (or Download the Files):**
    ```bash
    git clone https://github.com/linigoga/convert_Q3D_to_3D
    cd https://github.com/linigoga/convert_Q3D_to_3D
    ```
    If you don't use git, download the files (`transform_data.py`, `utils.pyx`, `setup.py`) and place them in the same directory.

2.  **Create a Virtual Environment (Recommended):**
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate  # On Linux/macOS
    # .venv\Scripts\activate   # On Windows
    conda activate my_env       # If using conda
    ```
    I recommend using using conda to create a virtual environment, as it is more reliable and easier to use.  If you don't have conda, you can install it with `pip install conda`.


3.  **Install Dependencies:**
    ```bash
    pip install numpy scipy h5py cython setuptools      # Using pip
    conda install numpy scipy h5py cython setuptools    # Using conda
    ```

4.  **Build the Cython Extension:**
    ```bash
    python setup.py build_ext --inplace
    ```
    This command compiles the `utils.pyx` file into a highly optimized C extension (or C++ extension, as configured) that can be imported by `transform_data.py`.

## 4. Usage

The tool is run from the command line using `transform_data.py`.

### Command-Line Arguments

`python transform_data.py INPUT DATA_TYPE [options]`

Convert all trimesteps fo electron charge densit, with all the modes:
`python transform_data.py /path/to/simulation/data charge -s electrons`

Convert all timesteps for the electric field, with all the modes:
`python transform_data.py /path/to/simulation/data field -f e2 -m 0 -t 100 200`

Convert Raw particle data from ions, for all timestpes specifying an output directory:
`python transform_data.py /path/to/simulation/data raw -s ions -o /path/to/output/directory`

Convert charge density for electrons, specifying a spatial domain, using a spatial average, and custom steps:
`python transform_data.py /path/to/data charge -s electrons -d 0.0 10.0 -5.0 5.0 -savg True -s_z 4 -s_r 4`

convert e1 field for a lineout, for modes 0 and 1:
`python transform_data.py /path/to/data field -f e1 -m 1 -line True`


## 10. Authors

- Lucas Ivan Iñigo Gamiz
- Bertrand Martinez
- Óscar Amaro
- Róbert Babjak