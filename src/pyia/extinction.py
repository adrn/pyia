# ruff: noqa

# Third-party
from __future__ import annotations

import numpy as np
import numpy.typing as npt

# https://github.com/mfouesneau/dustapprox/
# https://mfouesneau.github.io/dustapprox/


def get_ext_dr2_Babusiaux(
    G: npt.ArrayLike,
    bp: npt.ArrayLike,
    rp: npt.ArrayLike,
    ebv: npt.ArrayLike,
    maxnit: int = 8,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """Compute the Gaia extinctions assuming relations from Babusieux+2018
    Arguments: G, bp, rp, E(B-V)
    maxnit -- number of iterations
    Returns extinction in G, bp, rp
    Author: Sergey Koposov skoposov@cmu.edu
    """
    c1, c2, c3, c4, c5, c6, c7 = [
        0.9761,
        -0.1704,
        0.0086,
        0.0011,
        -0.0438,
        0.0013,
        0.0099,
    ]
    d1, d2, d3, d4, d5, d6, d7 = [
        1.1517,
        -0.0871,
        -0.0333,
        0.0173,
        -0.0230,
        0.0006,
        0.0043,
    ]
    e1, e2, e3, e4, e5, e6, e7 = [
        0.6104,
        -0.0170,
        -0.0026,
        -0.0017,
        -0.0078,
        0.00005,
        0.0006,
    ]
    A0 = 3.1 * ebv
    P1 = np.poly1d([c1, c2, c3, c4][::-1])

    def F1(bprp: npt.ArrayLike) -> npt.NDArray:
        return (
            np.poly1d([c1, c2, c3, c4][::-1])(bprp)
            + c5 * A0
            + c6 * A0**2
            + c7 * bprp * A0
        )

    def F2(bprp: npt.ArrayLike) -> npt.NDArray:
        return (
            np.poly1d([d1, d2, d3, d4][::-1])(bprp)
            + d5 * A0
            + d6 * A0**2
            + d7 * bprp * A0
        )

    def F3(bprp: npt.ArrayLike) -> npt.NDArray:
        return (
            np.poly1d([e1, e2, e3, e4][::-1])(bprp)
            + e5 * A0
            + e6 * A0**2
            + e7 * bprp * A0
        )

    xind = np.isfinite(bp + rp + G)
    curbp = bp - rp
    for i in range(maxnit):
        AG = F1(curbp) * A0
        Abp = F2(curbp) * A0
        Arp = F3(curbp) * A0
        curbp1 = bp - rp - Abp + Arp
        delta = np.abs(curbp1 - curbp)[xind]
        curbp = curbp1
    AG = F1(curbp) * A0
    Abp = F2(curbp) * A0
    Arp = F3(curbp) * A0
    return AG, Abp, Arp


__all__ = ["get_ext_dr2_Babusiaux"]
