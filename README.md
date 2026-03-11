<p align="center">
  <h1 align="center">sc-st-proAI Toolkit</h1>
  <p align="center"><strong>Discover reliable tools for single-cell, spatial-omics, and AI-driven biology</strong></p>
  <p align="center">
    <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License: MIT"></a>
    <img src="https://img.shields.io/badge/Python-3.12%2B-blue.svg" alt="Python 3.12+">
  </p>
</p>

A curated toolkit portal for single-cell transcriptomics, spatial-omics, multi-omics integration, and AI biology applications, deployed on GitHub Pages.

## 🌟 What is sc-st-proAI Toolkit?

In the rapidly evolving field of modern biology, finding reliable and well-maintained tools can be challenging. sc-st-proAI Toolkit provides a focused, data-driven leaderboard for key research areas, helping researchers avoid stagnant projects and find gold-standard solutions.

- **Focused Categories**: Single-cell, Spatial-omics, Multi-omics, ProteinAI, AI-Pathology
- **Purely Curated**: Domain-specific keyword matching with custom algorithms
- **Evidence-based**: Ranking based on popularity, growth, maintenance, and community impact

## Features

### 1. 📊 Tool Rankings
- Five specialized categories in modern biology
- Dual-track (Pipeline / Utility)
- Score breakdown showing popularity, trend, environment support, paper citations, and project health
- Project cards with preview images or initial avatars

### 2. 🤖 Automated Pipeline
- Weekly ranking refresh (scheduled via GitHub Actions)
- Automatic deployment to GitHub Pages

### 3. 💎 Data Features
- Install command detection
- Preview image extraction from README
- Project metadata and badge generation

### 4. 🎨 Frontend Portal
- Modern, responsive card design
- Easy category navigation
- Pipeline/Utility switching
- Mobile-friendly interface

## Project Structure

```
sc-st-proAI-toolkit/
├── .github/workflows/   # GitHub Actions automation
│   └── main.yml
├── scripts/             # Python crawler and ranking engine
│   └── sc_st_proai_gateway.py      # Core ranking engine
├── data/                # Data storage
│   ├── ranking_report.json
│   └── ranking_history.json
├── docs/                # GitHub Pages deployment
│   └── index.html
└── README.md
```

## Scoring Formulas

> ### Pipeline Score
> ```
> S = (Base + Trend + Env + Paper) × Health
> ```
> - **Base**: 6 × log10(Stars + 1)
> - **Trend**: 2.5 × log10(Weekly Growth + 1)
> - **Env**: +6 (Docker), +4 (Conda)
> - **Paper**: +4 if associated publication found
> - **Health**: Multiplier based on last update and issue health

> ### Utility Score
> ```
> S = (Base + Trend + Paper) × Health
> ```
> - **Base**: 7 × log10(Stars + 1)
> - **Trend**: 2.5 × log10(Weekly Growth + 1)
> - **Paper**: +4 if associated publication found
> - **Health**: Multiplier based on last update and issue health

## Quick Start

### Local Testing

```bash
# 1. Clone or download the repository
cd sc-st-proAI-toolkit

# 2. Install dependencies
pip install -r requirements.txt

# 3. Generate rankings (optional)
cd scripts
python sc_st_proai_gateway.py

# 4. Launch the website locally
cd ../docs
python -m http.server 8000
# Visit http://localhost:8000 in your browser
```

### Deployment

This project is designed to be deployed via GitHub Pages with automatic weekly updates:

1. Fork or use this repository as a template
2. Enable GitHub Pages in repository settings
3. The weekly GitHub Action will automatically update rankings

## Categories

- **Single-cell**: scRNA-seq analysis tools
- **Spatial-omics**: Spatial transcriptomics and visualization tools
- **Multi-omics**: Integration and multi-modal analysis tools
- **ProteinAI**: Protein structure prediction and design tools
- **AI-Pathology**: Computational pathology and digital pathology tools

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

Thanks to the broader open-source bioinformatics community for developing and maintaining these tools.

---

<p align="center">
  <i>Find reliable tools. Do great biology.</i>
</p>
