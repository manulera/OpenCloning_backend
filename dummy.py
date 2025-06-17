# from pydna.parsers import parse
from opencloning.assembly2 import Assembly, assembly2str, assemble
from pydna.dseqrecord import Dseqrecord


template = Dseqrecord('cccgaggggaatcgaa')
insert = Dseqrecord('ggggaatcAcccgag')

asm = Assembly((template, insert), limit=5, use_all_fragments=True)

# print(*asm.G.edges, sep='\n')


for a in asm.get_insertion_assemblies():
    if a[0][0] == 1:
        print(assembly2str(a))
        prod = assemble(asm.fragments, a, is_insertion=True)
        print(prod.seq)
        print()
