"""Exact sparse parent dilation; all current support/label semantics unchanged.

This is a separate qualification candidate, not the active production helper.
The sparse scatter is precisely the closed 3x3 Chebyshev neighborhood restricted
to the finite grid. Clipping an out-of-grid offset duplicates a boundary node
already reached by an in-grid offset; it cannot add an extra in-grid neighbor.
Every valid neighbor is reached by its original, unclipped offset.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from local_memory_logic import (
    EIGHT, overlap_incident_labels, support_components,
)


def parent_dilation_exact(previous: np.ndarray) -> np.ndarray:
    """Same values as binary_dilation(previous, EIGHT, border_value=0).

    Dense fallback is an implementation work choice, never a physical/support
    threshold: both branches compute the identical Boolean operation. Keep the
    original dense primitive when nine offsets per true parent exceed grid size.
    """
    previous = np.asarray(previous, dtype=bool)
    if previous.ndim != 2 or min(previous.shape) < 1:
        raise ValueError('parent mask must be a nonempty 2D grid')
    rows, cols = np.nonzero(previous)
    if 9 * len(rows) > previous.size:
        return ndimage.binary_dilation(previous, structure=EIGHT, border_value=0)
    result = np.zeros(previous.shape, dtype=bool)
    if len(rows) == 0:
        return result
    ny, nx = previous.shape
    row_offsets = [np.clip(rows+dy, 0, ny-1) for dy in (-1,0,1)]
    col_offsets = [np.clip(cols+dx, 0, nx-1) for dx in (-1,0,1)]
    for r in row_offsets:
        for c in col_offsets:
            result[r,c] = True
    return result


def next_suppression(phi: np.ndarray, previous_suppressed: np.ndarray | None,
                     radius: int, *, evaluate_overlap: bool,
                     include_component_sizes: bool = False) -> tuple[np.ndarray, dict]:
    """Original v002 assembly with only its parent dilation replaced.

    Original support_components and overlap detector are imported unchanged.
    Label ordering, inherited/incident IDs, component sizes and every output
    metadata value have the same meanings and calculation as v002.
    """
    labels, count, support = support_components(phi)
    if previous_suppressed is None:
        previous = np.zeros_like(support)
    else:
        previous = np.asarray(previous_suppressed, dtype=bool)
        if previous.shape != support.shape:
            raise ValueError('previous suppression field shape mismatch')
    parent_neighborhood = parent_dilation_exact(previous)
    inherited_values = np.unique(labels[parent_neighborhood & support])
    inherited_ids = [int(v) for v in inherited_values if int(v) != 0]
    if evaluate_overlap:
        incident = overlap_incident_labels(support, radius)
        overlap_centers = incident['overlap_centers']
        incident_ids = incident['incident_component_ids']
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
        'component_count_8': int(count),
        'component_sizes': sizes,
        'overlap_evaluated': bool(evaluate_overlap),
        'overlap_centers': overlap_centers,
        'inherited_component_ids': inherited_ids,
        'incident_component_ids': incident_ids,
        'suppressed_component_ids': suppressed_ids,
        'positive_P_nodes': int(np.count_nonzero(support)),
        'suppressed_support_nodes': int(np.count_nonzero(current_suppressed)),
    }
