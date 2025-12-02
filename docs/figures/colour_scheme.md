# rUTIs Colour Scheme

**Canonical Infection Status Colours**

This scheme is to be used for all plots where episodes or isolates are coloured by clinical infection status.

| Status | Hex Code | Colour Name |
| :--- | :--- | :--- |
| **UTI** | `#D55E00` | Vermilion (Warm Orange/Red) |
| **ASB** | `#0072B2` | Blue (Cool) |
| **Negative** | `#909090` | Grey (Neutral) |
| **Unknown** | `#CCCCCC` | Light Grey |

**Implementation in R**

Use the shared helpers in `R/plot_helpers.R`:

```r
source("R/plot_helpers.R")

ggplot(df, aes(x, y, color = Infection_Status)) +
  geom_point() +
  scale_colour_infection()
```

**Other Palettes**

*   **Culture-positive:** `#009E73` (Green)
*   **Within-Host:** `#0072B2`
*   **Between-Host:** `#CC79A7`
