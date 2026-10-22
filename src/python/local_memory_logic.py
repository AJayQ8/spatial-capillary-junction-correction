"""Versioned local support-overlap and Euler ancestry helpers.

The Boolean masks and label semantics match local_memory_logic.py. This version
avoids per-component full-grid scans in the per-Euler hot path.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

EIGHT = np.ones((3, 3), dtype=np.uint8)


def strict_product(phi: np.ndarray) -> np.ndarray:
    p = 27.0 * phi[0]
    p = p * phi[1]
    p = p * phi[2]
    return p


def support_components(phi: np.ndarray) -> tuple[np.ndarray, int, np.ndarray]:
    support = strict_product(phi) > 0.0
    labels, count = ndimage.label(support, structure=EIGHT)
    return labels, int(count), support


def overlap_incident_labels(support: np.ndarray, radius: int) -> dict:
    if support.dtype != np.dtype(bool) or support.ndim != 2:
        raise TypeError("support must be a 2D bool array")
    labels, count = ndimage.label(support, structure=EIGHT)
    if count < 2:
        return {"labels": labels, "count": int(count), "overlap_centers": np.zeros_like(support),
                "incident_component_ids": [], "overlap": False}
    width = 2 * int(radius) + 1
    nonzero = labels.copy()
    nonzero[nonzero == 0] = int(count) + 1
    overlap_centers = (ndimage.maximum_filter(labels, size=(width, width), mode="nearest") >
                       ndimage.minimum_filter(nonzero, size=(width, width), mode="nearest"))
    incident_neighborhood = ndimage.maximum_filter(overlap_centers.astype(np.uint8),
        size=(width, width), mode="nearest") > 0
    ids = np.unique(labels[support & incident_neighborhood])
    ids = [int(v) for v in ids if int(v) != 0]
    return {"labels": labels, "count": int(count), "overlap_centers": overlap_centers,
            "incident_component_ids": ids, "overlap": bool(ids)}


def next_suppression(phi: np.ndarray, previous_suppressed: np.ndarray | None,
                     radius: int, *, evaluate_overlap: bool,
                     include_component_sizes: bool = False) -> tuple[np.ndarray, dict]:
    labels, count, support = support_components(phi)
    if previous_suppressed is None:
        previous = np.zeros_like(support)
    else:
        previous = np.asarray(previous_suppressed, dtype=bool)
        if previous.shape != support.shape:
            raise ValueError("previous suppression field shape mismatch")
    parent_neighborhood = ndimage.binary_dilation(previous, structure=EIGHT, border_value=0)
    inherited_values = np.unique(labels[parent_neighborhood & support])
    inherited_ids = [int(v) for v in inherited_values if int(v) != 0]
    if evaluate_overlap:
        incident = overlap_incident_labels(support, radius)
        overlap_centers = incident["overlap_centers"]
        incident_ids = incident["incident_component_ids"]
    else:
        overlap_centers = np.zeros_like(support)
        incident_ids = []
    suppressed_ids = sorted(set(inherited_ids) | set(incident_ids))
    component_flags = np.zeros(count + 1, dtype=bool)
    if suppressed_ids:
        component_flags[np.asarray(suppressed_ids, dtype=np.intp)] = True
    current_suppressed = component_flags[labels] & support
    sizes = {}
    if include_component_sizes and count:
        counts = np.bincount(labels.ravel(), minlength=count + 1)
        sizes = {str(i): int(counts[i]) for i in range(1, count + 1)}
    return current_suppressed, {
        "component_count_8": int(count),
        "component_sizes": sizes,
        "overlap_evaluated": bool(evaluate_overlap),
        "overlap_centers": overlap_centers,
        "inherited_component_ids": inherited_ids,
        "incident_component_ids": incident_ids,
        "suppressed_component_ids": suppressed_ids,
        "positive_P_nodes": int(np.count_nonzero(support)),
        "suppressed_support_nodes": int(np.count_nonzero(current_suppressed)),
    }


def direct_overlap_oracle(support: np.ndarray, radius: int) -> tuple[bool, list[int]]:
    labels, count = ndimage.label(np.asarray(support, dtype=bool), structure=EIGHT)
    dilated = []
    for label_id in range(1, int(count) + 1):
        one = labels == label_id
        dilated.append(ndimage.maximum_filter(one.astype(np.uint8),
            size=(2 * radius + 1, 2 * radius + 1), mode="constant", cval=0) > 0)
    incident = set()
    for i in range(int(count)):
        for j in range(i + 1, int(count)):
            if np.any(dilated[i] & dilated[j]):
                incident.add(i + 1)
                incident.add(j + 1)
    return bool(incident), sorted(incident)

