import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar


@dataclass
class TempFileHandler:
    TempFiles: ClassVar[list] = []

    @staticmethod
    def write_temp_file(content: str) -> Path:
        with tempfile.NamedTemporaryFile(
            delete_on_close=False,
            delete=False,
            encoding="utf-8",
            suffix=".txt",
            mode="w+",
        ) as new_file:
            new_file.write(content)
            new_file.flush()
            TempFileHandler.TempFiles.append(new_file.name)
            return Path(new_file.name)

    @staticmethod
    def cleanup() -> None:
        for file in TempFileHandler.TempFiles:
            try:
                Path.unlink(file)
            except FileNotFoundError as e:
                print(f"Couldn't find temp file {file} for deletion: {e}")
