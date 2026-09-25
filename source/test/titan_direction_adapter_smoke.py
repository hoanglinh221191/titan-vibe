#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from titan_direction_adapter import (  # noqa: E402
    DirectionAdapterError,
    main,
    parse_vector_file,
)
from titan_path_controller import parse_key_value_file  # noqa: E402


class TitanDirectionAdapterSmoke(unittest.TestCase):
    def test_fixed_status(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            status = Path(temporary) / "status.txt"
            code = main(
                [
                    "--mode",
                    "fixed",
                    "--provider-id",
                    "straight-hj-a",
                    "--time-ps",
                    "20",
                    "--max-turn-angle-deg",
                    "20",
                    "--direction",
                    "1",
                    "2",
                    "3",
                    "--status",
                    str(status),
                ]
            )
            self.assertEqual(code, 0)
            values = parse_key_value_file(status)
            self.assertEqual(values["SCHEMA"], "TITAN_DIRECTION_STATUS_V1")
            self.assertEqual(values["DIRECTION_KIND"], "fixed")
            self.assertEqual(values["PROVIDER_ID"], "straight-hj-a")
            self.assertEqual(values["NEXT_VECTOR_Z_A"], "3")
            runtime = json.loads(
                (status.parent / "direction_provider_runtime.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                runtime["schema"], "titan-path-direction-provider-runtime-v1"
            )

    def test_vector_parser_rejects_zero_and_accepts_titan_line(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "DE_out_vector.dat"
            path.write_text(
                "0.2 0.3 0.4 Angle: 5 -- Modul: 0.54 Angstrom\n",
                encoding="utf-8",
            )
            self.assertEqual(parse_vector_file(path), (0.2, 0.3, 0.4))
            path.write_text("0 0 0\n", encoding="utf-8")
            with self.assertRaises(DirectionAdapterError):
                parse_vector_file(path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
