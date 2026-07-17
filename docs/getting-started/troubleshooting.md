# Troubleshooting

Common issues and solutions for ScanExitron users.

## Installation Issues

### Python Version Errors

??? question "ModuleNotFoundError or ImportError: No module named 'scanexitron'"

   **Symptom**: `ModuleNotFoundError: No module named 'scanexitron'`

   **Cause**: Wrong Python environment or installation failed.

   **Solution**:
   ```bash
   # Check Python version (must be >= 3.10 and < 3.13)
   python --version

   # Verify pip is using the correct Python
   which pip
   python -m pip --version

   # Reinstall in the current environment
   python -m pip install --force-reinstall scanexitron
   ```

### Missing External Tools

??? question "FileNotFoundError or subprocess error related to regtools/samtools/bedtools"

   **Symptom**: `Required tool 'regtools' not found. Please install it and retry.`

   **Cause**: The external bioinformatics binaries are missing or not in your PATH.

   **Solution**:
   - If using Conda: Re-run `conda install -c bioconda scanexitron` to ensure dependencies are resolved correctly.
   - If installing via pip/source: Make sure the paths to `regtools`, `samtools`, and `bedtools` are added to your environment `PATH`:
     ```bash
     export PATH=$PATH:/path/to/binaries/
     ```

## Runtime Issues

### MAPQ Filtering Failures

??? question "MAPQ filtering failed. Exiting."

   **Symptom**: ScanExitron stops early during the MAPQ filter step.

   **Cause**: BAM file may not be indexed, or index file is missing, or BAM is corrupted.

   **Solution**:
   Ensure the BAM file index (`.bai`) is present in the **same directory** as the `.bam` file.
   ```bash
   # If index is missing, generate it:
   samtools index input.bam
   ```

## General Help

### Enable Verbose Logging

For debugging, enable detailed logging output:

```bash
scanexitron run ... --verbose
```

### Check Version

Ensure you're using the latest version:

```bash
# Check current version
scanexitron --version
```

## Getting Further Help

If your issue isn't covered here:

1. **Check existing issues**: [GitHub Issues](https://github.com/ylab-hi/ScanExitron/issues)
2. **Open a new issue**: Include:
   - ScanExitron version (`scanexitron --version`)
   - Python version (`python --version`)
   - Operating system
   - Complete error message/traceback
   - Minimal reproducible example/description of the steps taken
