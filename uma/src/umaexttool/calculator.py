from __future__ import annotations

import warnings

with warnings.catch_warnings():
    # suppress the missing PySisiphus warning
    warnings.simplefilter("ignore")
    from fairchem.core import pretrained_mlip, FAIRChemCalculator


def init(model: str, device: str) -> FAIRChemCalculator:
    """Initialize an UmaCalculator instance."""
    predictor=pretrained_mlip.get_predict_unit("uma-s-1", device=device)
    return FAIRChemCalculator(predictor, task_name=model)
