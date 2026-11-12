# Zenodo / GitHub release checklist

The repository is packaged in `06_repository/`. To obtain the DOI cited in
the manuscript:

1. Create the GitHub repository
   ```
   cd 06_repository
   git init
   # copy in: 00_pipeline/  02_scripts/  03_results/tables/supplementary/
   #          04_manuscript/supplementary/  Snakefile  environment.yml
   #          Dockerfile  RUN.md  README.md  LICENSE  CITATION.cff
   git add -A && git commit -m "v1.0.0 — manuscript release"
   git remote add origin https://github.com/luwen001199/neurospatial-gradients.git
   git push -u origin main
   ```
   Recommended `.gitignore`:
   ```
   01_data/
   03_results/figures/*.png
   03_results/tables/*.csv
   *.pyc
   __pycache__/
   ```
   (Figures and tables regenerate; ship only the supplementary tables that
   the manuscript cites.)

2. Tag the release
   ```
   git tag -a v1.0.0 -m "First release accompanying the manuscript"
   git push origin v1.0.0
   ```

3. Connect the repository to Zenodo (zenodo.org → GitHub → enable the repo),
   then open a release on GitHub. Zenodo mints the DOI automatically and
   deposits the tagged snapshot.

4. Fill in the placeholders before releasing:
   - `[OWNER]` in `README.md` and `CITATION.cff`
   - author names and ORCIDs in `CITATION.cff`
   - the assigned DOI in the manuscript's Data & Code Availability section
     (currently `10.5281/zenodo.[TO-BE-ASSIGNED]`).

5. After the DOI exists, update `04_manuscript/manuscript_main.md`
   (Section 7) and `CITATION.cff`, tag `v1.0.1`, and let Zenodo version the
   DOI.
