"""
PNR discrimination histograms for a P-SNSPD readout.

Computes, per 16-sample waveform, six discrimination observables:
  1. max_amp        - peak amplitude
  2. area           - area under curve (trapezoid)
  3. area_sq        - area under squared curve (trapezoid)   <- energy-like
  4. slew           - max single-sample rising slope         <- Cahall-style
  5. max_amp_sinc   - peak amplitude after sinc interpolation
  6. area_sq_sinc   - area under squared curve after sinc interpolation

Processing is chunked so the full ~50M-waveform dataset never has to be
upsampled in one allocation.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import resample, resample_poly

# np.trapz (NumPy <2) was renamed np.trapezoid (NumPy >=2)
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

# ----------------------------------------------------------------------------
# Config -- set these to match your acquisition
# ----------------------------------------------------------------------------
WF_LEN        = 16          # samples per waveform
UPSAMPLE      = 16          # sinc interpolation factor (16 -> 256 pts/wf)
N_BASE        = 4           # leading samples used for baseline estimate (0 = off)
POLARITY      = +1          # +1 if pulses go up, -1 if they go down
OFFSET_BINARY = False       # True if ADC data is offset binary (XOR 0x8000)
SINC_METHOD   = "fft"       # "fft" (resample) or "poly" (resample_poly)
CHUNK         = 1_000_000   # waveforms processed per batch


# ----------------------------------------------------------------------------
# Core processing
# ----------------------------------------------------------------------------
def _to_float(chunk_int):
    """Offset-binary -> signed if needed, then float64."""
    if OFFSET_BINARY:
        chunk_int = (chunk_int.astype(np.uint16) ^ 0x8000).astype(np.int16)
    return chunk_int.astype(np.float64)


def _baseline_correct(w):
    """Subtract per-waveform baseline (mean of first N_BASE samples) and
    apply polarity so pulses are positive-going."""
    if N_BASE > 0:
        w = w - w[:, :N_BASE].mean(axis=1, keepdims=True)
    return POLARITY * w


def _sinc_upsample(w):
    """Band-limited interpolation onto a WF_LEN*UPSAMPLE grid.

    fft:  exact Whittaker-Shannon on a periodic window (fast, but assumes the
          window is one period -> some ringing if the pulse hasn't decayed by
          the last sample).
    poly: polyphase windowed-sinc; more robust to non-periodic edges.
    """
    n_out = WF_LEN * UPSAMPLE
    if SINC_METHOD == "fft":
        return resample(w, n_out, axis=1)
    return resample_poly(w, UPSAMPLE, 1, axis=1)


def compute_observables(wfs, chunk=CHUNK, verbose=True):
    """wfs: (N, WF_LEN) array of raw samples. Returns dict of (N,) float arrays."""
    N = wfs.shape[0]
    out = {k: np.empty(N, dtype=np.float64) for k in
           ("max_amp", "area", "area_sq", "slew", "max_amp_sinc", "area_sq_sinc")}

    for s in range(0, N, chunk):
        e = min(s + chunk, N)
        w = _baseline_correct(_to_float(wfs[s:e]))

        out["max_amp"][s:e] = w.max(axis=1)
        out["area"][s:e]    = _trapz(w, axis=1)
        out["area_sq"][s:e] = _trapz(w ** 2, axis=1)
        out["slew"][s:e]    = np.diff(w, axis=1).max(axis=1)

        up = _sinc_upsample(w)
        out["max_amp_sinc"][s:e] = up.max(axis=1)
        # dx = 1/UPSAMPLE keeps the integral on the same physical scale as area_sq
        out["area_sq_sinc"][s:e] = _trapz(up ** 2, axis=1, dx=1.0 / UPSAMPLE)

        del w, up
        if verbose:
            print(f"  processed {e:,}/{N:,}", end="\r")
    if verbose:
        print()
    return out


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
def _hist(ax, vals, title, bins=2000):
    """Histogram with log bins when the data is strictly positive and spans a
    decent dynamic range; linear bins otherwise (e.g. signed slew / area)."""
    vals = vals[np.isfinite(vals)]
    use_log = vals.min() > 0 and (vals.max() / vals.min()) > 50
    if use_log:
        edges = np.logspace(np.log10(vals.min()), np.log10(vals.max()), bins)
    else:
        edges = np.linspace(vals.min(), vals.max(), bins)

    ax.hist(vals, bins=edges, histtype="step", color="C0", lw=1.0)
    if use_log:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.set_ylabel("counts")


def plot_observables(obs, suptitle="PNR discrimination", savepath=None):
    panels = [
        ("max_amp",      "Max amplitude"),
        ("area",         "Area under curve"),
        ("area_sq",      "Area under curve (squared)"),
        ("slew",         "Slew rate (max rising slope)"),
        ("max_amp_sinc", "Max amplitude (sinc interp.)"),
        ("area_sq_sinc", "Area sq. (sinc interp.)"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, (key, title) in zip(axes.ravel(), panels):
        _hist(ax, obs[key], title)
    fig.suptitle(suptitle)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=150, bbox_inches="tight")
    return fig


# ----------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    data = np.load('/Users/gabriel/received_data_RAM_0.5.npy')   # your data
    wfs = data.reshape(-1, WF_LEN)[:1_000_000]
    print(f"{wfs.shape[0]:,} waveforms of length {WF_LEN}")

    obs = compute_observables(wfs)
    plot_observables(obs, savepath="pnr_discrimination.png")
    plt.show()
