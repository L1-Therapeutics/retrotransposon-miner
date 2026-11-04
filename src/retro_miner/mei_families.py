"""Family tokens shared by the caller, RepeatMasker filter, and plots.

Human calls stay ALU / LINE1 / SVA. Mouse adds the active classes that
RepeatMasker does not file under those three names: B1, B2, IAP, and ETn.
B1 is checked before Alu because mouse B1 consensuses are ``SINE/Alu``.
"""

from __future__ import annotations

import re

MEI_CALL_FAMILIES = frozenset({"ALU", "LINE1", "SVA", "B1", "B2", "IAP", "ETN"})

_B1_TOKEN = re.compile(r"(?:^|[^A-Z0-9])P?B1(?:F[0-9]*|D[0-9]+)?(?:[^A-Z0-9]|$)")
_B2_TOKEN = re.compile(r"(?:^|[^A-Z0-9])B(?:2|3A?)(?:[^A-Z0-9]|$)")
# ETn is a prefix glued to the rest of the name (ETnERV3), so the
# trailing boundary used for B1/B2 would miss it.
_ETN_TOKEN = re.compile(r"(?:^|[^A-Z0-9])ETN")


def normalize_mei_family_token(token: str) -> str:
    """Map a consensus or RepeatMasker label onto a call family, or ``""``."""
    t = (token or "").upper()
    if not t:
        return ""
    if _B1_TOKEN.search(t) or "SINE/B1" in t:
        return "B1"
    if _B2_TOKEN.search(t) or "SINE/B2" in t:
        return "B2"
    if "IAP" in t:
        return "IAP"
    if _ETN_TOKEN.search(t):
        return "ETN"
    if "SVA" in t:
        return "SVA"
    if "ALU" in t:
        return "ALU"
    if "LINE1" in t or "L1" in t:
        return "LINE1"
    return ""
