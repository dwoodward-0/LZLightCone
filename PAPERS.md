# Third-party documents used by this analysis

These five public PDFs are declared inputs in `astra.yaml`, under `data/`. Their
content is **not** distributed with this repository: the files appear as
git-annex links with no content behind them. To run the analysis, download each
one to the path shown and check it against the SHA256 below. The digitization
scripts verify these hashes and stop if a file doesn't match.

| File (`data/…`) | ASTRA input | Document | Download | SHA256 |
|---|---|---|---|---|
| `lz_ws2024_paper.pdf` | `lz_ws2024_paper` | LZ Collaboration, *Dark Matter Search Results from 4.2 Tonne-Years of Exposure of the LUX-ZEPLIN (LZ) Experiment*, arXiv:2410.17036v3 | <https://arxiv.org/pdf/2410.17036v3> | `d358bffbaed7eaddd431f9b540790673d650f00d0964282a089a15c84b2f9544` |
| `lz_ws2024_tevpa_talk.pdf` | `lz_ws2024_tevpa_talk` | LZ Collaboration, WS2024 results talk, TeVPA 2024 | <https://indico.uchicago.edu/event/427/contributions/1325/attachments/359/548/lz_results_tevpa.pdf> | `617bddb94385af704021fb495c676ea81fb6103bcddf51b392843b80e27301f6` |
| `lz_backgrounds_paper.pdf` | `lz_backgrounds_paper` | LZ Collaboration, *Background Determination for the LUX-ZEPLIN (LZ) Dark Matter Experiment*, arXiv:2211.17120 | <https://arxiv.org/abs/2211.17120> | `52223054fcdcb65c14787ec4c82c5a2685d222cfd189aa29bdbc230829c7c2e7` |
| `lz_er_efficiency_paper.pdf` | `lz_er_efficiency_paper` | LZ Collaboration, public ER efficiencies, arXiv:2511.17350v1 | <https://arxiv.org/pdf/2511.17350v1> | `9c4b10d7c93dcbb4323d49c9a892b6b21e9e4c3d7ba55d14eeba6064330d2738` |
| `lz_sensitivity2018_paper.pdf` | `lz_sensitivity2018_paper` | LZ Collaboration, *Projected WIMP Sensitivity of the LUX-ZEPLIN (LZ) Dark Matter Experiment*, arXiv:1802.06039v2 | <https://arxiv.org/pdf/1802.06039v2> | `5913d97d2d45710e4a6c7fddf8aac3a98f20adc3be3faec7c553a98611f9c058` |

Check a download with `shasum -a 256 data/<file>.pdf`. If the arXiv version
links stop matching, an earlier or later version was fetched.

The literature quoted as evidence in `astra.yaml` (the original NEST paper,
the NEST review, the statistics-conventions paper and others) is cached by
`astra-tools` outside the repository. Fetch it with
`uvx astra-tools@0.2.18 paper` as needed.

Other files in `data/` come from HEPData records (cited in `astra.yaml`) or
were derived from the documents above by the `src/digitize_*.py` scripts. They
are uploaded with the repository.
