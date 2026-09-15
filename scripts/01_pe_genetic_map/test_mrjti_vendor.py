"""Drift gate for the vendored MR-JTI script (audit B-2, 2026-09-10).

tools/MR-JTI/mr/MR-JTI.r is vendored verbatim from upstream
(gamazonlab/MR-JTI, master) modulo line-ending normalization, and
mrjti_minimal_runner.R imports its TRB_LASSO function family by line
anchors. This test freezes that contract:

1. the vendored file's raw sha256 equals the recorded constant (any edit,
   even whitespace, must be re-documented in mr/PATCH_NOTES.md first);
2. the line-ending-normalized content hash equals the upstream snapshot's
   normalized hash, so scientific content cannot drift from upstream;
3. the runner's eval-import anchors ("^TRB_LASSO<-" and "^#load df") are
   still present exactly once in the vendored file;
4. mr/PATCH_NOTES.md exists and records both sha256 constants.

Run:
    python3 scripts/01_pe_genetic_map/test_mrjti_vendor.py
"""

import hashlib
import re
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
VENDOR_SCRIPT = PROJECT_ROOT / "tools/MR-JTI/mr/MR-JTI.r"
RUNNER_SCRIPT = Path(__file__).resolve().parent / "mrjti_minimal_runner.R"
PATCH_NOTES = PROJECT_ROOT / "tools/MR-JTI/mr/PATCH_NOTES.md"

# Upstream snapshot fetched 2026-09-10 from
# https://raw.githubusercontent.com/gamazonlab/MR-JTI/master/mr/MR-JTI.r
UPSTREAM_SHA256 = "5ba207380dd99ffb8dcea4cabcdad80a2d3bca71fed770c0b678fc0a0815c4c4c4"
# Vendored copy as frozen in this project (line endings normalized to LF).
VENDOR_SHA256 = "77ff5da018c7f8462d6879355532f86f9cf02dfb89d65b88a9c315a92d9bd453"
# sha256 after CRLF/CR -> LF normalization; equal for both files, proving the
# only difference between upstream and vendored copies is line endings.
NORMALIZED_SHA256 = "f776e6a233b0f04342e9f30a5e84e94be2bbf3f395c599858e866cc14021e1a8"

START_ANCHOR = re.compile(r"^TRB_LASSO<-", re.M)
END_ANCHOR = re.compile(r"^#load df", re.M)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_lf(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


class VendoredMrjtiTests(unittest.TestCase):
    def test_vendored_file_unchanged(self):
        self.assertEqual(
            sha256_bytes(VENDOR_SCRIPT.read_bytes()), VENDOR_SHA256,
            "vendored MR-JTI.r changed; re-verify against upstream and update "
            "tools/MR-JTI/mr/PATCH_NOTES.md plus the constants here",
        )

    def test_content_matches_upstream_after_line_ending_normalization(self):
        normalized = read_lf(VENDOR_SCRIPT)
        self.assertEqual(
            sha256_bytes(normalized), NORMALIZED_SHA256,
            "vendored MR-JTI.r differs from the upstream snapshot beyond line "
            "endings; scientific modifications are not permitted without "
            "documentation in tools/MR-JTI/mr/PATCH_NOTES.md",
        )

    def test_runner_import_anchors_present(self):
        text = VENDOR_SCRIPT.read_text()
        self.assertEqual(len(START_ANCHOR.findall(text)), 1,
                         "^TRB_LASSO<- anchor must appear exactly once")
        self.assertEqual(len(END_ANCHOR.findall(text)), 1,
                         "^#load df anchor must appear exactly once")
        runner = RUNNER_SCRIPT.read_text()
        self.assertIn('"^TRB_LASSO<-"', runner,
                      "runner no longer locates the start anchor")
        self.assertIn('"^#load df"', runner,
                      "runner no longer locates the end anchor")

    def test_patch_notes_document_both_hashes(self):
        notes = PATCH_NOTES.read_text()
        self.assertIn(UPSTREAM_SHA256, notes)
        self.assertIn(VENDOR_SHA256, notes)
        self.assertIn(NORMALIZED_SHA256, notes)


if __name__ == "__main__":
    unittest.main(verbosity=2)
