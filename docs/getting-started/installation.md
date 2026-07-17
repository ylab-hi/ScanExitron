# Installation

Get ScanExitron installed on your system.

## Prerequisites

!!! info "Requirements"
    - **Python**: 3.10 to 3.12 (Python 3.13 is currently not supported due to compatibility restrictions in external dependencies)
    - **Operating System**: Linux, macOS
    - **External Tools**: The following external tools must be installed and available on your system `PATH` (if not installing via Conda):
        - **regtools**: `0.5.0`
        - **samtools**: `≥ 1.10`
        - **bedtools**: `≥ 2.26`

## Installation Methods

Choose your preferred installation method:

=== "Bioconda (Recommended)"

    The easiest way to install ScanExitron, as it automatically installs all required external command-line dependencies (like `regtools`, `samtools`, and `bedtools`):

    ```bash
    conda install -c bioconda scanexitron
    ```

=== "pip"

    Install the Python package from PyPI. Note that you must install the external dependencies (`regtools`, `samtools`, and `bedtools`) separately.

    ```bash
    pip install scanexitron
    ```

    Verify the installation:

    ```bash
    scanexitron --version
    ```

=== "uv"

    Using the fast `uv` package manager:

    ```bash
    uv pip install scanexitron
    ```

    Or create a virtual environment first:

    ```bash
    uv venv --python 3.11
    source .venv/bin/activate
    uv pip install scanexitron
    ```

=== "from source"

    For development or to build the latest features from GitHub:

    ```bash
    # Clone the repository
    git clone https://github.com/ylab-hi/ScanExitron.git
    cd ScanExitron

    # Install package and development dependencies
    pip install -e ".[dev]"

    # Verify installation
    scanexitron --version
    ```

## Verification

Confirm ScanExitron is installed correctly:

```bash
# Check version
scanexitron --version

# View available commands
scanexitron --help
```

## Troubleshooting Installation

### Common Issues

??? question "ImportError: No module named 'scanexitron'"

    **Solution**: Ensure you have activated the correct Python environment:

    ```bash
    # Check which Python is being used
    which python

    # Reinstall in the current environment
    pip install --force-reinstall scanexitron
    ```

??? question "External dependency 'regtools' not found"

    **Solution**: Ensure `regtools` is installed and on your system `PATH`. If you didn't install via Conda, you can download and build regtools from source:

    ```bash
    git clone https://github.com/griffithlab/regtools.git
    cd regtools
    mkdir build
    cd build
    cmake ..
    make
    # Add the build directory to your PATH
    export PATH=$PATH:$(pwd)
    ```

For more issues, see the [Troubleshooting Guide](troubleshooting.md).

## Next Steps

Now that ScanExitron is installed, try the [Quick Start](quick-start.md) tutorial to run your first calling job!
