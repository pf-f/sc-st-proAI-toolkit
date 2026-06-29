# sc-st-proAI Toolkit

Discover reliable open-source tools for single-cell analysis, spatial omics, multi-omics integration, protein AI, and AI pathology.

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![Python 3.12+](https://img.shields.io/badge/Python-3.12%2B-blue.svg)

## Live Site

Open the deployed toolkit here:

**https://pf-f.github.io/sc-st-proAI-toolkit/**

## At A Glance

- Static GitHub Pages portal in `docs/`
- Python crawler and ranking engine in `scripts/`
- Weekly automated ranking refresh through GitHub Actions
- Public JSON outputs in `data/` and `docs/data/`

## What It Does

sc-st-proAI Toolkit is a GitHub Pages portal backed by an automated crawler and ranking engine. It collects public GitHub projects, filters for biology and biomedical AI relevance, classifies projects into focused categories, and presents ranked Pipeline and Utility leaderboards.

The current public site is intended for:

- finding maintained workflows and libraries quickly
- comparing environment support such as Docker and Conda
- spotting paper-linked and actively updated projects
- discovering installation commands and README preview images

## Categories

- Single-cell
- Spatial-omics
- Multi-omics
- ProteinAI
- AI-Pathology

## Ranking Signals

Pipeline scores consider:

- GitHub stars
- weekly star growth
- Docker and Conda support
- paper or citation evidence
- repository freshness and issue health

Utility scores use similar signals, with stronger emphasis on project popularity and maintenance.

## Crawler Improvements

The crawler searches across complementary GitHub query styles:

- exact keyword matches in repository name, description, and README
- topic-based searches for category aliases
- recent activity and minimum-star thresholds

It then reads README content before final filtering, so projects are not rejected only because their description or topics are sparse. Category matching is centralized to avoid drift between the collector and report generator.

Runtime knobs:

```bash
SEARCH_CUTOFF_DATE=2025-01-01
MIN_STARS=10
SEARCH_MAX_PAGES=5
REQUEST_TIMEOUT=30
WEBSITE_URL=https://pf-f.github.io/sc-st-proAI-toolkit/
```

## Local Development

```bash
pip install -r requirements.txt

cd scripts
python sc_st_proai_gateway.py

cd ../docs
python -m http.server 8000
```

Then open `http://localhost:8000`.

## Deployment

The repository is configured for GitHub Pages from the `docs/` folder. The deployed URL is:

```text
https://pf-f.github.io/sc-st-proAI-toolkit/
```

The weekly GitHub Action refreshes rankings, copies generated JSON into `docs/data/`, commits changed data, and lets GitHub Pages publish the static site.

To trigger a refresh manually, run the `sc-st-proAI Toolkit Weekly Update` workflow from the GitHub Actions tab.

## Project Structure

```text
sc-st-proAI-toolkit/
  .github/workflows/main.yml
  data/
    ranking_report.json
    ranking_history.json
  docs/
    index.html
    data/
      ranking_report.json
      ranking_history.json
  scripts/
    sc_st_proai_gateway.py
  requirements.txt
```

## License

MIT License. See [LICENSE](LICENSE).
