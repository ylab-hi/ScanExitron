# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.4.0] - 2026-07-16

### Added
- Expanded test coverage for `core`, `gtf`, and `vcf` modules.
- GitHub Actions workflow for linting, type-checking, and running tests across Python 3.10-3.12.
- Bioconda recipe (`conda/meta.yaml`).

### Changed
- Merged `exitron2vcf` into `scanexitron convert`, using an `--input` flag.
- Limited supported Python versions to `>=3.10,<3.13` due to a `regtools`/zlib constraint.
- Eliminated `config.ini` in favor of CLI arguments.

[1.4.0]: https://github.com/ylab-hi/ScanExitron/releases/tag/v1.4.0
