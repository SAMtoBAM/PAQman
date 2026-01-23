# Telomeric sequence strings/regex to use for different species

These basic strings have been confirmed by us and users to detect telomeres in assemblies/chromosomes manually verified to contains telomeric ends.
PLEASE FEEL FREE TO PUSH ADDITIONS IF YOU CONFIRM THE DETECTION

## Yeast
_Saccharomyces cerevisiae_ = "GGTGTG" or "TG{1,3}TG{1,3}" <br/>
_Schizosaccharomyces japonicus_ = "T{0,1}A{2,4}C{2,4}T{0,1}A{1,2}GAC" or "G[TG][CG]T{2,3}A" <br/>
_Candida albicans_ = "ACTTCTTGGTGTACGGATGTCTA"

## Filamentous fungi
_Fusarium_ = "TTAGGG" (default) <br/>
_Magnaportha/Pyricularia_ = "TTAGGG" (default) <br/>

## Animals
Human = "TTAGGG" (default)

## Plants
_Oryza sativa_ (rice) = "TTTAGGG"

