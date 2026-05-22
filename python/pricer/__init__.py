"""pricer — option pricing & risk engine with a C++ core.

Thin Python package over the compiled ``_pricer`` extension. Example::

    import pricer
    c = pricer.black_scholes_call(100, 100, 0.05, 0.20, 1.0)
    g = pricer.black_scholes_greeks(pricer.OptionType.Call, 100, 100, 0.05, 0.20, 1.0)
    iv = pricer.implied_vol(pricer.OptionType.Call, c, 100, 100, 0.05, 1.0)
"""
from ._pricer import (  # noqa: F401
    OptionType,
    Greeks,
    RiskMeasures,
    black_scholes_price,
    black_scholes_call,
    black_scholes_put,
    black_scholes_greeks,
    implied_vol,
    mc_price,
    mc_price_parallel,
    qmc_price,
    var_es,
    __version__,
)

__all__ = [
    "OptionType",
    "Greeks",
    "RiskMeasures",
    "black_scholes_price",
    "black_scholes_call",
    "black_scholes_put",
    "black_scholes_greeks",
    "implied_vol",
    "mc_price",
    "mc_price_parallel",
    "qmc_price",
    "var_es",
    "__version__",
]
