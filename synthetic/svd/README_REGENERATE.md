# Regenerating the layered synthetic SVD dataset

The generated arrays are not shipped (they total more than 1 GB). From this directory:

```bash
python make_synthetic_svd_dataset.py
python make_synthetic_svd_layered_tensor.py
```

See `README_synthetic_svd_levels.md` and `README_layered_tensor.md` for the intended layout and outputs.
