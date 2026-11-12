import sys
sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
import run_sensitivity as RS
out, store = RS.run_metabric_robustness(n_boot=500, seed=42)
print("\nDONE")
