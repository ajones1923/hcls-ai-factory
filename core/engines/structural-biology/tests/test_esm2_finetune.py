"""B11: ESM-2 LoRA fine-tuning. Real run (loads esm2_t12_35M, trains on a synthetic signal)."""
import random
import pytest
from esm2_finetune import ESM2FineTuner

AAS = "ACDEFGHIKLMNPQRSTVWY"


@pytest.mark.slow
def test_lora_finetune_learns_and_is_efficient():
    random.seed(0)
    seqs = ["".join(random.choice(AAS) for _ in range(40)) for _ in range(24)]
    labels = [float(s.count("K")) for s in seqs]          # learnable: lysine count
    ft = ESM2FineTuner(task="regression")
    res = ft.finetune(seqs, labels, epochs=10, lr=5e-4)
    assert res["improved"]                                 # loss decreased
    assert res["trainable_pct"] < 10.0                     # LoRA: only a few % trainable
    preds = ft.predict(seqs[:3])
    assert len(preds) == 3 and all(isinstance(p, float) for p in preds)
