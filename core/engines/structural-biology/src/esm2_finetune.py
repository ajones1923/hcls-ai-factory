"""
ESM-2 domain fine-tuning via LoRA (B11).

Parameter-efficient fine-tuning of a protein language model for a custom property
(regression or classification): a small ESM-2 backbone + a LoRA adapter + a task head, trained
on the GB10 GPU. LoRA keeps the trainable parameter count tiny (a few %), so domain models train
fast and ship small — the path to bespoke property predictors (solubility, activity, stability)
without the external developability repos.
"""
from __future__ import annotations

from typing import Any


class ESM2FineTuner:
    def __init__(self, model_name: str = "facebook/esm2_t12_35M_UR50D",
                 task: str = "regression", num_labels: int = 1) -> None:
        self.model_name = model_name
        self.task = task
        self.num_labels = num_labels
        self._tok = None
        self._model = None
        self._torch = None
        self.device = "cpu"

    def _build(self) -> None:
        if self._model is not None:
            return
        import torch
        from transformers import AutoTokenizer, EsmForSequenceClassification
        from peft import LoraConfig, get_peft_model, TaskType
        self._torch = torch
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self._tok = AutoTokenizer.from_pretrained(self.model_name)
        problem = "regression" if self.task == "regression" else "single_label_classification"
        base = EsmForSequenceClassification.from_pretrained(
            self.model_name, num_labels=self.num_labels, problem_type=problem)
        lora = LoraConfig(task_type=TaskType.SEQ_CLS, r=8, lora_alpha=16,
                          target_modules=["query", "value"], lora_dropout=0.05)
        self._model = get_peft_model(base, lora).to(self.device)

    def param_stats(self) -> dict[str, int]:
        self._build()
        total = sum(p.numel() for p in self._model.parameters())
        trainable = sum(p.numel() for p in self._model.parameters() if p.requires_grad)
        return {"trainable": trainable, "total": total,
                "trainable_pct": round(100.0 * trainable / total, 3)}

    def finetune(self, sequences: list[str], labels: list[float], epochs: int = 3,
                 lr: float = 1e-3) -> dict[str, Any]:
        self._build()
        torch = self._torch
        opt = torch.optim.AdamW((p for p in self._model.parameters() if p.requires_grad), lr=lr)
        enc = self._tok(sequences, return_tensors="pt", padding=True, truncation=True,
                        max_length=512).to(self.device)
        y = torch.tensor(labels, dtype=torch.float32, device=self.device).unsqueeze(1) \
            if self.task == "regression" else \
            torch.tensor(labels, dtype=torch.long, device=self.device)
        self._model.train()
        loss_curve = []
        for _ in range(epochs):
            opt.zero_grad()
            out = self._model(**enc, labels=y)
            out.loss.backward()
            opt.step()
            loss_curve.append(round(float(out.loss.item()), 5))
        return {"epochs": epochs, "device": self.device, "loss_curve": loss_curve,
                "improved": loss_curve[-1] < loss_curve[0], **self.param_stats()}

    def predict(self, sequences: list[str]) -> list[float]:
        self._build()
        torch = self._torch
        self._model.eval()
        enc = self._tok(sequences, return_tensors="pt", padding=True, truncation=True,
                        max_length=512).to(self.device)
        with torch.no_grad():
            logits = self._model(**enc).logits
        return logits.squeeze(-1).cpu().tolist() if self.task == "regression" \
            else logits.argmax(-1).cpu().tolist()
