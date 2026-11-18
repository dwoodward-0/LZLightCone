"""Decision flags for the statistics recipes (PLR, FC, their figures).

Each recipe passes its declared decisions as CLI flags; a script accepts only the
options it implements and refuses the rest, so a universe asking for an
unimplemented option fails loudly instead of silently running the default.
"""
from __future__ import annotations

IMPLEMENTED = {
    "statistics": {
        "plr": ("both_with_crosscheck", "profile_likelihood_2d"),
        "fc": ("both_with_crosscheck", "cut_and_count_fc"),
        "both": ("both_with_crosscheck",),
    },
    "background_components": ("all_table1",),
    # combined_single = one pooled likelihood sample; the realized exposure is the
    # WS2024 3.3 t-y of the Table-I counts and signal archive (see AGENTS.md).
    "exposure_structure": ("combined_single",),
    "median_sensitivity_definition": ("median_of_toy_limits",),
    "power_constraint": ("applied_1sigma",),
    "test_statistic_sidedness": ("two_sided",),
}


def add_decision_args(ap, names, method):
    """Add --<decision> flags (required) restricted to the implemented options."""
    for name in names:
        choices = IMPLEMENTED[name][method] if name == "statistics" else IMPLEMENTED[name]
        ap.add_argument("--" + name.replace("_", "-"), required=True, choices=choices,
                        help=f"decision {name} (implemented: {', '.join(choices)})")


def decisions_of(args, names):
    return {n: getattr(args, n) for n in names}
