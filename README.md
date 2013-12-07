metapynets
==========

Metabolic pathways aligment using Petri nets. `metapynets` is an implementation  of algorithm described on the
paper [Local piecewise alignment of metabolic pathways][1]. It allows compare two metabolic pathways represented as Petri
nets in PNML (Petri Net Markup Language).

- Install dependencies: 

`pip install -r requirements.txt`

- Usage:

`python metapynets.py pnml/pnml_hsa00010.xml pnml/pnml_mmu00010.xml`

[1]: http://bulma.net/~marctc/PosterBiotecno.pdf "Local piecewise alignment of metabolic pathways"
