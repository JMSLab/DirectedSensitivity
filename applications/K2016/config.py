import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent

with open(HERE / 'config.toml', 'rb') as _stream:
    CONFIG = tomllib.load(_stream)

PATHS = {name: str(HERE / path) for name, path in CONFIG['paths'].items()}


def output(name):
    return str(Path(PATHS['output_dir']) / name)
