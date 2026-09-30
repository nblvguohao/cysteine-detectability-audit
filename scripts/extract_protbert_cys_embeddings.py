"""Extract frozen ProtBERT embeddings for 31-aa Cys-centred site windows on CUDA."""
from __future__ import annotations

import argparse
import csv
import json
import platform
import re
import sys
import time
from pathlib import Path

import numpy as np
import torch
import transformers
from transformers import AutoModel, AutoTokenizer

from common import INPUTS, ROOT, load_frozen, sha256, write_json


OUT = ROOT / "plm"
DEFAULT_MODEL = ROOT / "models" / "prot_bert"
VALID = set("ACDEFGHIKLMNPQRSTVWXY")


def load_sequences(path: Path):
    sequences = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sequences[row["Entry"]] = row["Sequence"].strip().upper()
    return sequences


def normalize(sequence: str):
    return "".join(residue if residue in VALID else "X" for residue in sequence)


def site_window(sequence: str, one_based_position: int, flank: int = 15):
    index = one_based_position - 1
    if index < 0 or index >= len(sequence) or sequence[index] != "C":
        raise ValueError(f"Expected C at position {one_based_position}, length={len(sequence)}")
    start = max(0, index - flank)
    stop = min(len(sequence), index + flank + 1)
    residues = normalize(sequence[start:stop])
    center_residue_index = index - start
    return " ".join(residues), center_residue_index + 1  # account for [CLS]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--batch-size", type=int, default=128)
    args = parser.parse_args()
    started = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    X, y, folds, proteins, components, positions, _ = load_frozen()
    sequence_path = INPUTS / "sly_proteome.tsv"
    sequences = load_sequences(sequence_path)
    missing = sorted(set(proteins) - set(sequences))
    if missing:
        raise ValueError(f"Missing {len(missing)} benchmark proteins, examples={missing[:10]}")

    windows, center_tokens = [], []
    normalized_residue_count = 0
    for accession, position in zip(proteins, positions):
        raw = sequences[accession]
        window, center = site_window(raw, int(position))
        normalized_residue_count += sum(char not in VALID for char in raw[max(0, int(position)-16):int(position)+15])
        windows.append(window)
        center_tokens.append(center)

    tokenizer = AutoTokenizer.from_pretrained(
        args.model_dir, local_files_only=True, do_lower_case=False
    )
    model = AutoModel.from_pretrained(args.model_dir, local_files_only=True)
    model.eval().to("cuda:0")
    torch.cuda.reset_peak_memory_stats()
    embeddings = np.empty((len(windows), int(model.config.hidden_size)), dtype=np.float16)

    with torch.inference_mode():
        for start in range(0, len(windows), args.batch_size):
            stop = min(len(windows), start + args.batch_size)
            encoded = tokenizer(
                windows[start:stop],
                add_special_tokens=True,
                padding=True,
                return_tensors="pt",
            )
            encoded = {key: value.to("cuda:0") for key, value in encoded.items()}
            indices = torch.as_tensor(center_tokens[start:stop], device="cuda:0")
            with torch.autocast(device_type="cuda", dtype=torch.float16):
                hidden = model(**encoded).last_hidden_state
            selected = hidden[torch.arange(stop - start, device="cuda:0"), indices]
            embeddings[start:stop] = selected.float().cpu().numpy().astype(np.float16)
            if start == 0:
                center_ids = encoded["input_ids"][torch.arange(stop - start, device="cuda:0"), indices]
                center_symbols = tokenizer.convert_ids_to_tokens(center_ids.cpu().tolist())
                if set(center_symbols) != {"C"}:
                    raise ValueError(f"Tokenizer/site alignment failed: {set(center_symbols)}")
            if stop % (args.batch_size * 25) == 0 or stop == len(windows):
                print(json.dumps({"embedded_sites": stop, "total_sites": len(windows)}), flush=True)

    output = OUT / "protbert_cys_embeddings.npz"
    np.savez(
        output,
        embeddings=embeddings,
        proteins=proteins,
        positions=positions,
        folds=folds,
        y=y,
    )
    write_json(OUT / "protbert_embedding_audit.json", {
        "completed": True,
        "representation": "frozen ProtBERT last-layer Cys token from variable-length 31-aa centred window",
        "flank_residues": 15,
        "boundary_handling": "shorten window at termini; no artificial padding residues",
        "ambiguous_residue_handling": "map residues outside ACDEFGHIKLMNPQRSTVWXY to X",
        "normalized_residue_occurrences_in_windows": int(normalized_residue_count),
        "n_sites": int(len(windows)),
        "n_proteins": int(len(set(proteins))),
        "embedding_shape": list(embeddings.shape),
        "embedding_dtype": str(embeddings.dtype),
        "batch_size": int(args.batch_size),
        "elapsed_seconds": time.time() - started,
        "peak_gpu_memory_bytes": int(torch.cuda.max_memory_allocated()),
        "hardware": {
            "platform": platform.platform(),
            "gpu": torch.cuda.get_device_name(0),
            "capability": list(torch.cuda.get_device_capability(0)),
        },
        "versions": {
            "python": sys.version,
            "numpy": np.__version__,
            "torch": torch.__version__,
            "torch_cuda": torch.version.cuda,
            "transformers": transformers.__version__,
        },
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                INPUTS / "frozen_benchmark_data.npz",
                INPUTS / "primary_site_folds.csv",
                sequence_path,
                args.model_dir / "config.json",
                args.model_dir / "vocab.txt",
                args.model_dir / "pytorch_model.bin",
            ]
        },
        "output_sha256": sha256(output),
    })


if __name__ == "__main__":
    main()
