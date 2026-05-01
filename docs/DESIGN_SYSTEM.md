# Design System: Publication Visual Standards

## 🎨 Color Palette
We use a curated, high-contrast palette to distinguish between different categories of genetic evidence.

| Category | Hex Code | Usage |
| :--- | :--- | :--- |
| **Novel Locus** | `#d4a62a` | Gold border/point to highlight new discoveries. |
| **Known Locus** | `#243046` | Dark slate for established genetic associations. |
| **High Confidence** | `#215ca8` | Deep blue for targets with LAVA/COLOC support. |
| **Immune-Mediated** | `#d98a29` | Orange for PheWAS traits in the immune category. |
| **Non-Immune** | `#8b8f97` | Grey for secondary pleiotropic associations. |

## 📊 Typography scale
- **Main Titles**: 21pt - 24pt (DejaVu Serif, Bold)
- **Subtitles**: 11pt - 12pt (DejaVu Sans, Regular)
- **Axis Labels**: 10pt - 11pt (DejaVu Sans)
- **Point Annotations**: 9pt (DejaVu Sans)

## 📐 Layout Rules
- **Manhattan Plots**: Vertical scale should be $-log_{10}(P)$ with a red dashed line at $5 \times 10^{-8}$.
- **Pleiotropy Maps**: Linear horizontal layout to prevent text overlap in multi-trait figures.
- **Forest Plots**: Odds Ratio (OR) on X-axis, Log-scale preferred for highly skewed effects.
- **Supplementary Figures**: Must use the "Minimalist B&W" style (White boxes, Black text, no background colors).

## ✅ Do's and Don'ts
| Do | Don't |
| :--- | :--- |
| Use `cowplot` for multi-panel figures. | Use generic red/green colors for status. |
| Embed fonts in PDFs for publication. | Leave overlapping labels in locus plots. |
| Use high-DPI (300+) PNG for presentations. | Use "METAL" as a software name in legends. |
| Include sample sizes in figure legends. | Use default R colors (red1, blue1). |

## 🖼️ Component reference
- **Lollipop Plots**: Used for ranking loci by significance.
- **Lane Plots**: Used for therapeutic mapping (Gene -> Target -> Drug).
- **Dot Plots**: Used for pathway enrichment bubbles (Size = Count, Color = Significance).
