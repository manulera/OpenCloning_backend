from pydna.parsers import parse
from opencloning.assembly2 import Assembly, assembly2str
from pydna.dseqrecord import Dseqrecord

template = parse('NC_000913.gb')[0]
insert = parse('rph_homology_arm.gb')[0]

template = Dseqrecord('cccgaggggaatcgaa')
insert = Dseqrecord('Acccgagggggaatc')


asm = Assembly((template, insert), limit=5, use_all_fragments=True)

# print(*asm.G.edges, sep='\n')

for a in asm.get_insertion_assemblies():
    print(assembly2str(a))


# hom = 'CAATTTG'

# template = Dseqrecord('cccgaggggaat')
# insert1 = Dseqrecord(f'cgaggggTTAA{hom}')
# insert2 = Dseqrecord(f'{hom}CCAggggaa')


# asm = Assembly((template, insert1, insert2), limit=5, use_all_fragments=True)

# print(*asm.G.edges, sep='\n')

# for a in asm.get_insertion_assemblies():
#     print(assembly2str(a))
# for a in asm.assemble_insertion():
#     print(a.seq)
