from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from typing import Any

from RNAmotiFold.Misc.EnumTypes import SubOptMode
from RNAmotiFold.Misc.Helpers import interval_check


@dataclass(frozen=True)
class SubOptEnergy:
    subopt_energy: pos_decimal | None = None
    subopt_mode: SubOptMode = SubOptMode.PERCENT

    @property
    def commandline_parameter(self):
        match self.subopt_mode:
            case SubOptMode.ABSOLUTE:
                return "-e"
            case SubOptMode.PERCENT:
                return "-c"


class limited_temperature(float):
    def __new__(cls, value: Any) -> float:
        val = float(value)
        if interval_check(0, 100, val):
            return super().__new__(cls, val)
        raise ValueError("Temperature can only be between 0°C and 100°C")


class fraction(float):
    def __new__(cls, value: Any) -> float:
        val = float(value)
        if interval_check(0, 1, val):
            return super().__new__(cls, val)
        raise ValueError("Probability filter can only be between 0 and 1")


class pos_int(int):
    def __new__(cls, value: Any) -> int:
        val = int(value)
        if val > 0:
            return super().__new__(cls, val)
        raise ValueError("Only positive integers allowed here!")


class pos_float(float):
    def __new__(cls, value: Any) -> float:
        val = float(value)
        if val >= 0:
            return super().__new__(cls, val)
        raise ValueError("Only positive floats allowed here!")


class pos_decimal(Decimal):
    def __new__(cls, value: Any) -> Decimal:
        val = Decimal(value)
        if val >= 0:
            return super().__new__(cls, val)
        raise ValueError("Only positive decimals allowed here!")
