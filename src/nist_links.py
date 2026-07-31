from __future__ import annotations

from urllib.parse import quote_plus

from .matching import normalize_cas


def nist_gc_url(cas_number: object, chemical_name: object, mask: int = 2000) -> str:
    cas = normalize_cas(cas_number)
    if cas:
        return f"https://webbook.nist.gov/cgi/cbook.cgi?ID=C{cas}&Mask={mask}"
    return f"https://webbook.nist.gov/cgi/cbook.cgi?Name={quote_plus(str(chemical_name or ''))}&Units=SI&Mask={mask}"
