# Telomeric sequence strings/regex to use for different species

These basic strings have been confirmed by us and users to detect telomeres in assemblies/chromosomes manually verified to contains telomeric ends. <br/>
Please push additions/open issues if you have any confirmed repeats and share the data so I can give it whirl!

| Lineage/group  | species/genus |  repeat | alternate repeat |
| -------------- | ------- | ------- | ---------------- |
| Yeast | _Saccharomyces cerevisiae_ | GGTGTG | TG{1,3}TG{1,3} |
| Yeast | _Schizosaccharomyces japonicus_ | T{0,1}A{2,4}C{2,4}T{0,1}A{1,2}GAC | G[TG][CG]T{2,3}A |
| Yeast | _Candida albicans_ | ACTTCTTGGTGTACGGATGTCTA | TG{1,3}TG{1,3} |
|  |  |  |  |
| Filamentous fungi | _Fusarium_ | TTAGGG |  |
| Filamentous fungi | _Magnaportha/Pyricularia_ | TTAGGG |  |
|  |  |  |  |
| Animal | Human | TTAGGG |  |
|  |  |  |  |
| Plants | _Oryza sativa_ | TTTAGGG |  |
