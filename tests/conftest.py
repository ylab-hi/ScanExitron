from pathlib import Path
import pytest

TEST_DATA_DIR = Path(__file__).parent


@pytest.fixture
def test_bam() -> Path:
    return TEST_DATA_DIR / "test.bam"


@pytest.fixture
def test_exitron() -> Path:
    return TEST_DATA_DIR / "test.exitron"
