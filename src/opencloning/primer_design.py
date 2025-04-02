from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design, assembly_fragments
from Bio.SeqFeature import SimpleLocation
from pydna.utils import locations_overlap, shift_location, location_boundaries
from pydna.amplicon import Amplicon
from .pydantic_models import PrimerModel
from Bio.Seq import reverse_complement
from Bio.Restriction.Restriction import RestrictionType
from Bio.Data.IUPACData import ambiguous_dna_values as _ambiguous_dna_values
from primer3.bindings import calc_tm as _calc_tm
from typing import Callable


def primer3_calc_tm(seq: str) -> float:
    return _calc_tm(seq.upper())


ambiguous_dna_values = _ambiguous_dna_values.copy()
# Remove acgt
for base in 'ACGT':
    del ambiguous_dna_values[base]


def homologous_recombination_primers(
    pcr_seq: Dseqrecord,
    pcr_loc: SimpleLocation,
    hr_seq: Dseqrecord,
    hr_loc: SimpleLocation,
    homology_length: int,
    minimal_hybridization_length: int,
    insert_forward: bool,
    target_tm: float,
    spacers: list[str] | None = None,
    tm_func: Callable[[str], float] = primer3_calc_tm,
    estimate_function: Callable[[str], float] | None = None,
) -> tuple[str, str]:

    fragment2amplify = pcr_loc.extract(pcr_seq)
    amplicon = primer_design(
        fragment2amplify,
        limit=minimal_hybridization_length,
        target_tm=target_tm,
        tm_func=tm_func,
        estimate_function=estimate_function,
    )

    if insert_forward:
        fwd_primer, rvs_primer = amplicon.primers()
    else:
        rvs_primer, fwd_primer = amplicon.primers()

    if fwd_primer is None or rvs_primer is None:
        raise ValueError('Primers could not be designed, try changing settings.')

    hr_loc_start, hr_loc_end = location_boundaries(hr_loc)
    fwd_homology_start = hr_loc_start - homology_length
    rvs_homology_end = hr_loc_end + homology_length

    if spacers is None:
        spacers = ['', '']

    if len(spacers) != 2:
        raise ValueError("The 'spacers' list must contain exactly two elements.")

    if not hr_seq.circular:
        if fwd_homology_start < 0:
            raise ValueError('Forward homology region is out of bounds.')
        if rvs_homology_end > len(hr_seq):
            raise ValueError('Reverse homology region is out of bounds.')

    # Convert to locations
    fwd_arm = shift_location(SimpleLocation(fwd_homology_start, hr_loc_start), 0, len(hr_seq))
    rvs_arm = shift_location(SimpleLocation(hr_loc_end, rvs_homology_end), 0, len(hr_seq))

    if locations_overlap(fwd_arm, rvs_arm, len(hr_seq)):
        raise ValueError('Homology arms overlap.')

    fwd_homology = fwd_arm.extract(hr_seq)
    rvs_homology = rvs_arm.extract(hr_seq).reverse_complement()

    fwd_primer_seq = (str(fwd_homology.seq) + spacers[0]).lower() + str(fwd_primer.seq).upper()
    rvs_primer_seq = (str(rvs_homology.seq) + reverse_complement(spacers[1])).lower() + str(rvs_primer.seq).upper()

    return fwd_primer_seq, rvs_primer_seq


def gibson_assembly_primers(
    templates: list[Dseqrecord],
    homology_length: int,
    minimal_hybridization_length: int,
    target_tm: float,
    circular: bool,
    spacers: list[str] | None = None,
    tm_func: Callable[[str], float] = primer3_calc_tm,
    estimate_function: Callable[[str], float] | None = None,
) -> list[PrimerModel]:

    initial_amplicons = [
        primer_design(
            template,
            limit=minimal_hybridization_length,
            target_tm=target_tm,
            tm_func=tm_func,
            estimate_function=estimate_function,
        )
        for template in templates
    ]

    for i, amplicon in enumerate(initial_amplicons):
        if None in amplicon.primers():
            raise ValueError(f'Primers could not be designed for template {templates[i].name}, try changing settings.')

    maxlink = 40
    if spacers is not None:
        maxlink = max(len(spacer) for spacer in spacers)
        spacers = [Dseqrecord(spacer) for spacer in spacers]
        initial_amplicons_with_spacers = []
        # For linear assemblies, the first spacer is the first thing
        if not circular:
            initial_amplicons_with_spacers.append(spacers.pop(0))
        for amplicon in initial_amplicons:
            initial_amplicons_with_spacers.append(amplicon)
            initial_amplicons_with_spacers.append(spacers.pop(0))
        initial_amplicons = initial_amplicons_with_spacers
        # Maxlink is used to define what is a spacer or what is a template (see docs)

    assembly_amplicons: list[Amplicon] = assembly_fragments(
        initial_amplicons, overlap=homology_length, circular=circular, maxlink=maxlink
    )

    all_primers = sum((list(amplicon.primers()) for amplicon in assembly_amplicons), [])

    for i in range(0, len(all_primers), 2):
        fwd, rvs = all_primers[i : i + 2]
        template = templates[i // 2]
        template_name = template.name if template.name != 'name' else f'seq_{template.id}'
        fwd.name = f'{template_name}_fwd'
        rvs.name = f'{template_name}_rvs'

    return [PrimerModel(id=0, name=primer.name, sequence=str(primer.seq)) for primer in all_primers]


def sanitize_enzyme_site(site: str) -> str:
    """
    Replaces any ambiguous bases with the first value of the ambiguous base from the Biopython IUPACData.
    """
    return ''.join(b if b not in ambiguous_dna_values else ambiguous_dna_values[b][0] for b in site)


def simple_pair_primers(
    template: Dseqrecord,
    minimal_hybridization_length: int,
    target_tm: float,
    left_enzyme: RestrictionType,
    right_enzyme: RestrictionType,
    filler_bases: str,
    spacers: list[str] | None = None,
    left_enzyme_inverted: bool = False,
    right_enzyme_inverted: bool = False,
    tm_func: Callable[[str], float] = primer3_calc_tm,
    estimate_function: Callable[[str], float] | None = None,
) -> tuple[PrimerModel, PrimerModel]:
    """
    Design primers to amplify a DNA fragment, if left_enzyme or right_enzyme are set, the primers will be designed
    to include the restriction enzyme sites.
    """

    if spacers is None:
        spacers = ['', '']

    if len(spacers) != 2:
        raise ValueError("The 'spacers' list must contain exactly two elements.")

    amplicon = primer_design(
        template,
        limit=minimal_hybridization_length,
        target_tm=target_tm,
        tm_func=tm_func,
        estimate_function=estimate_function,
    )
    fwd_primer, rvs_primer = amplicon.primers()

    if fwd_primer is None or rvs_primer is None:
        raise ValueError('Primers could not be designed, try changing settings.')

    template_name = template.name if template.name != 'name' else f'seq_{template.id}'

    left_site = '' if left_enzyme is None else sanitize_enzyme_site(left_enzyme.site)
    right_site = '' if right_enzyme is None else sanitize_enzyme_site(right_enzyme.site)

    if left_enzyme_inverted:
        left_site = reverse_complement(left_site)
    if right_enzyme_inverted:
        right_site = reverse_complement(right_site)
    if left_enzyme is not None:
        left_site = filler_bases + left_site
    if right_enzyme is not None:
        right_site = filler_bases + right_site

    fwd_primer_seq = left_site + spacers[0] + fwd_primer.seq
    rvs_primer_seq = right_site + reverse_complement(spacers[1]) + rvs_primer.seq

    fwd_primer_name = f'{template_name}_{left_enzyme}_fwd' if left_enzyme is not None else f'{template_name}_fwd'
    rvs_primer_name = f'{template_name}_{right_enzyme}_rvs' if right_enzyme is not None else f'{template_name}_rvs'

    return (
        PrimerModel(id=0, name=fwd_primer_name, sequence=str(fwd_primer_seq)),
        PrimerModel(id=0, name=rvs_primer_name, sequence=str(rvs_primer_seq)),
    )


# def gateway_attB_primers(
#     template: Dseqrecord,
#     minimal_hybridization_length: int,
#     target_tm: float,
#     sites: tuple[str, str],
#     spacers: tuple[str, str],
#     filler_bases: str = 'GGGG',
# ) -> tuple[PrimerModel, PrimerModel]:
#     if spacers is None:
#         spacers = ['', '']

#     if len(spacers) != 2:
#         raise ValueError("The 'spacers' list must contain exactly two elements.")

#     if sites[0] not in primer_design_attB or sites[1] not in primer_design_attB:
#         raise ValueError('Invalid attB site.')

#     amplicon = primer_design(template, limit=minimal_hybridization_length, target_tm=target_tm)
#     fwd_primer, rvs_primer = amplicon.primers()

#     if fwd_primer is None or rvs_primer is None:
#         raise ValueError('Primers could not be designed, try changing settings.')

#     template_name = template.name if template.name != 'name' else f'seq_{template.id}'

#     left_site = primer_design_attB[sites[0]]
#     right_site = primer_design_attB[sites[1]]

#     fwd_primer_seq = filler_bases + left_site + spacers[0] + fwd_primer.seq
#     rvs_primer_seq = filler_bases + right_site + reverse_complement(spacers[1]) + rvs_primer.seq

#     return (
#         PrimerModel(id=0, name=f'{template_name}_{sites[0]}_fwd', sequence=str(fwd_primer_seq)),
#         PrimerModel(id=0, name=f'{template_name}_{sites[1]}_rvs', sequence=str(rvs_primer_seq)),
#     )
