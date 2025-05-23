"""Slightly different assembly implementation"""

from pydna.utils import (
    shift_location as _shift_location,
    flatten,
    location_boundaries as _location_boundaries,
    locations_overlap as _locations_overlap,
)
from pydna._pretty import pretty_str as _pretty_str
from pydna.common_sub_strings import common_sub_strings as common_sub_strings_str
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna.dseq import Dseq as _Dseq
from pydna.primer import Primer as _Primer
from pydna.seqrecord import SeqRecord as _SeqRecord
import networkx as _nx
import itertools as _itertools
from Bio.SeqFeature import SimpleLocation, Location
from .dna_utils import sum_is_sticky
from Bio.Seq import reverse_complement
from Bio.Restriction.Restriction import RestrictionBatch, AbstractCut
import regex
import copy

# Currently unused, commented out because it's not tested
# def primers_clash(assembly, fragments):
#     edge_pairs = zip(assembly, assembly[1:])
#     for (_u1, _v1, _, start_location), (_u2, _v2, end_location, _) in edge_pairs:
#         # Only for primer joins
#         if not isinstance(fragments[abs(_v1) - 1], _Dseqrecord):
#             continue
#         if _locations_overlap(start_location, end_location, len(fragments[abs(_v1) - 1])):
#             return True
#     return False


def limit_iterator(iterator, limit):
    for i, x in enumerate(iterator):
        if i >= limit:
            raise ValueError(f'Too many possible paths (more than {limit})')
        yield x


def gather_overlapping_locations(locs: list[Location], fragment_length: int):
    """
    Turn a list of locations into a list of tuples of those locations, where each tuple contains
    locations that overlap. For example, if locs = [loc1, loc2, loc3], and loc1 and loc2 overlap,
    the output will be [(loc1, loc2), (loc3,)].
    """
    # Make a graph with all the locations as nodes
    G = _nx.Graph()
    for i, loc in enumerate(locs):
        G.add_node(i, location=loc)

    # Add edges between nodes that overlap
    for i in range(len(locs)):
        for j in range(i + 1, len(locs)):
            if _locations_overlap(locs[i], locs[j], fragment_length):
                G.add_edge(i, j)

    # Get groups of overlapping locations
    groups = list()
    for loc_set in _nx.connected_components(G):
        groups.append(tuple(locs[i] for i in loc_set))

    # Sort by location of the first element in each group (does not matter which since they are overlapping)
    groups.sort(key=lambda x: _location_boundaries(x[0])[0])

    return groups


# def assembly_checksum(G: _nx.MultiDiGraph, edge_list):
#     """Calculate a checksum for an assembly, from a list of edges in the form (u, v, key)."""
#     checksum_list = list()
#     for edge in edge_list:
#         u, v, key = edge
#         checksum_list.append(G.get_edge_data(u, v, key)['uid'])

#     return min('-'.join(checksum_list), '-'.join(checksum_list[::-1]))


def ends_from_cutsite(cutsite: tuple[tuple[int, int], AbstractCut], seq: _Dseq):
    if cutsite is None:
        raise ValueError('None is not supported')

    cut_watson, cut_crick, ovhg = seq.get_cut_parameters(cutsite, is_left=None)
    if ovhg < 0:
        # TODO check the edge in circular
        return (
            ("5'", str(seq[cut_watson:cut_crick].reverse_complement()).lower()),
            ("5'", str(seq[cut_watson:cut_crick]).lower()),
        )
    elif ovhg > 0:
        return (
            ("3'", str(seq[cut_crick:cut_watson]).lower()),
            ("3'", str(seq[cut_crick:cut_watson].reverse_complement()).lower()),
        )

    return ('blunt', ''), ('blunt', '')


def restriction_ligation_overlap(
    seqx: _Dseqrecord, seqy: _Dseqrecord, enzymes=RestrictionBatch, partial=False, allow_blunt=False
):
    """Find overlaps. Like in stiky and gibson, the order matters"""
    cuts_x = seqx.seq.get_cutsites(*enzymes)
    cuts_y = seqy.seq.get_cutsites(*enzymes)
    # If blunt ends are allowed, something similar to this could be done to allow
    # joining with linear sequence ends, but for now it messes up with the only_adjacent_edges
    # case
    # if allow_blunt:
    #     if not seqx.circular:
    #         cuts_x.append(((len(seqx), 0), None))
    #     if not seqy.circular:
    #         cuts_y.append(((0, 0), None))
    matches = list()
    for cut_x, cut_y in _itertools.product(cuts_x, cuts_y):
        # A blunt end
        if allow_blunt and cut_x[0][1] == cut_y[0][1] == 0:
            matches.append((cut_x[0][0], cut_y[0][0], 0))
            continue

        # Otherwise, test overhangs
        overlap = sum_is_sticky(ends_from_cutsite(cut_x, seqx.seq)[0], ends_from_cutsite(cut_y, seqy.seq)[1], partial)
        if not overlap:
            continue
        x_watson, x_crick, x_ovhg = seqx.seq.get_cut_parameters(cut_x, is_left=False)
        y_watson, y_crick, y_ovhg = seqy.seq.get_cut_parameters(cut_y, is_left=True)
        # Positions where the overlap would start for full overlap
        left_x = x_watson if x_ovhg < 0 else x_crick
        left_y = y_watson if y_ovhg < 0 else y_crick

        # Correct por partial overlaps
        left_x += abs(x_ovhg) - overlap

        matches.append((left_x, left_y, overlap))
    return matches


def combine_algorithms(*algorithms):
    """Combine algorithms, if any of them returns a match, the match is returned."""

    def combined(seqx, seqy, limit):
        matches = list()
        for algorithm in algorithms:
            matches += algorithm(seqx, seqy, limit)
        return matches

    return combined


def blunt_overlap(seqx: _Dseqrecord, seqy: _Dseqrecord, limit=None):
    """Find blunt overlaps"""
    if seqx.seq.three_prime_end()[0] == 'blunt' and seqy.seq.five_prime_end()[0] == 'blunt':
        return [(len(seqx), 0, 0)]
    return []


def common_sub_strings(seqx: _Dseqrecord, seqy: _Dseqrecord, limit=25):
    return common_sub_strings_str(str(seqx.seq).upper(), str(seqy.seq).upper(), limit)


def gibson_overlap(seqx: _Dseqrecord, seqy: _Dseqrecord, limit=25):
    """
    The order matters, we want alignments like:

    oooo------xxxx
              xxxx------oooo
    Product: oooo------xxxx------oooo

    Not like:

              oooo------xxxx
    xxxx------oooo
    Product (unwanted): oooo
    """

    # Because Gibson enzymes remove 5' overhangs, we remove them from the sequence
    # when looking for homology, then we shift the location of the second fragment accordingly.
    # This is only relevant for linear fragments, so we don't need to worry about
    # shifting locations for circular fragments.
    trim_x_left = -seqx.seq.ovhg if seqx.seq.ovhg < 0 else 0
    trim_x_right = seqx.seq.watson_ovhg() if seqx.seq.watson_ovhg() < 0 else None
    trim_y_left = -seqy.seq.ovhg if seqy.seq.ovhg < 0 else 0
    trim_y_right = seqy.seq.watson_ovhg() if seqy.seq.watson_ovhg() < 0 else None

    stringx = str(seqx.seq[trim_x_left:trim_x_right]).upper()
    stringy = str(seqy.seq[trim_y_left:trim_y_right]).upper()
    # We have to convert to list because we need to modify the matches
    matches = [
        list(m) for m in common_sub_strings_str(stringx, stringy, limit) if (m[1] == 0 and m[0] + m[2] == len(stringx))
    ]
    for match in matches:
        match[0] += trim_x_left
        match[1] += trim_y_left

    # convert to tuples again
    return [tuple(m) for m in matches]


def sticky_end_sub_strings(seqx: _Dseqrecord, seqy: _Dseqrecord, limit=0):
    """For now, if limit 0 / False only full overlaps are considered."""
    overlap = sum_is_sticky(seqx.seq.three_prime_end(), seqy.seq.five_prime_end(), limit)
    if overlap:
        return [(len(seqx) - overlap, 0, overlap)]
    return []


def zip_match_leftwards(seqx: _SeqRecord, seqy: _SeqRecord, match: tuple[int, int, int]):
    """Starting from the rightmost edge of the match, return a new match encompassing the max
    number of bases. This can be used to return a longer match if a primer aligns for longer
    than the limit or a shorter match if there are mismatches. This is convenient to maintain
    as many features as possible.

    >>> seq = _Dseqrecord('AAAAACGTCCCGT')
    >>> primer = _Dseqrecord('ACGTCCCGT')
    >>> match = (13, 9, 0) # an empty match at the end of each
    >>> zip_match_leftwards(seq, primer, match)
    (4, 0, 9)

    Works in circular molecules if the match spans the origin:
    >>> seq = _Dseqrecord('TCCCGTAAAAACG', circular=True)
    >>> primer = _Dseqrecord('ACGTCCCGT')
    >>> match = (6, 9, 0)
    >>> zip_match_leftwards(seq, primer, match)
    >>> (10, 0, 9)

    """

    query_x = seqrecord2str_for_alignment(seqx)
    query_y = seqrecord2str_for_alignment(seqy)

    # In circular sequences, the match may go beyond the left-most edge of the sequence if it spans
    # the origin:
    # Primer:          ACGTCCCGT
    #                  |||||||||
    # Circular seq:    ACGTCCCGT -> Equivalent to Dseqrecord('CCCGTACGT', circular=True)
    #                      ^
    #                      Origin
    # We would start from the last T and move leftwards, but we would stop at the origin
    # For those cases we shift by length, then go back

    end_on_x = match[0] + match[2]
    if isinstance(seqx, _Dseqrecord) and seqx.circular and end_on_x <= len(seqx):
        end_on_x += len(seqx)

    end_on_y = match[1] + match[2]
    if isinstance(seqy, _Dseqrecord) and seqy.circular and end_on_y <= len(seqy):
        end_on_y += len(seqy)

    count = 0
    for x, y in zip(reversed(query_x[:end_on_x]), reversed(query_y[:end_on_y])):
        if x != y:
            break
        count += 1

    # Shift back by length if needed
    start_on_x = (end_on_x - count) % len(seqx)
    start_on_y = (end_on_y - count) % len(seqy)

    return (start_on_x, start_on_y, count)


def zip_match_rightwards(seqx: _Dseqrecord, seqy: _Dseqrecord, match: tuple[int, int, int]):
    """Same as zip_match_leftwards, towards the right."""

    query_x = seqrecord2str_for_alignment(seqx)
    query_y = seqrecord2str_for_alignment(seqy)

    start_on_x, start_on_y, _ = match
    count = 0
    for x, y in zip(query_x[start_on_x:], query_y[start_on_y:]):
        if x != y:
            break
        count += 1
    return (start_on_x, start_on_y, count)


def seqrecord2str_for_alignment(seqr: _SeqRecord):
    """Transform a Dseqrecord to a string representation where U is replaced by T, everything is upper case and
    circular sequences are repeated twice."""
    out = str(seqr.seq).upper().replace('U', 'T')
    if isinstance(seqr, _Dseqrecord) and seqr.circular:
        return out * 2
    return out


def alignment_sub_strings(seqx: _Dseqrecord | _Primer, seqy: _Dseqrecord | _Primer, limit=25, mismatches=0):
    """"""

    if isinstance(seqx, _Primer) and isinstance(seqy, _Dseqrecord):
        primer = seqx
        template = seqy
        reverse_primer = False
    elif isinstance(seqx, _Dseqrecord) and isinstance(seqy, _Primer):
        primer = seqy
        template = seqx
        reverse_primer = True
    else:
        raise ValueError('One of the sequences must be a primer and the other a Dseqrecord')

    if len(primer) < limit:
        return []

    subject = seqrecord2str_for_alignment(template)
    query = (
        seqrecord2str_for_alignment(primer[:limit]) if reverse_primer else seqrecord2str_for_alignment(primer[-limit:])
    )

    re_matches = list(regex.finditer('(' + query + '){s<=' + str(mismatches) + '}', subject, overlapped=True))
    re_matches += list(regex.finditer('(?r)(' + query + '){s<=' + str(mismatches) + '}', subject, overlapped=True))

    out = set()
    for re_match in re_matches:

        start, end = re_match.span()

        # For circular sequences the same match is returned twice unless it falls
        # on the origin, we eliminate duplicates here
        if start >= len(template):
            continue

        # This extends match beyond the limit if the primer aligns more than that
        # and reduces the match if the primer has mismatches
        if reverse_primer:
            # Match in the same format as other assembly algorithms
            starting_match = (start, 0, end - start)
            out.add(zip_match_rightwards(template, primer, starting_match))
        else:
            # Match in the same format as other assembly algorithms
            starting_match = (len(primer) - limit, start, end - start)
            out.add(zip_match_leftwards(primer, template, starting_match))

    return list(sorted(out))


def fill_left(seq: _Dseq):
    """Fill the left overhang of a sequence with the complementary sequence."""
    new_watson = seq.watson
    new_crick = seq.crick

    # Watson 5' overhang
    if seq.ovhg < 0:
        new_crick = new_crick + reverse_complement(seq.watson[: -seq.ovhg])
    # Crick 5' overhang
    elif seq.ovhg > 0:
        new_watson = reverse_complement(seq.crick[-seq.ovhg :]) + new_watson

    return _Dseq(new_watson, new_crick, 0)


def fill_right(seq: _Dseq):
    """Fill the right overhang of a sequence with the complementary sequence."""
    new_watson = seq.watson
    new_crick = seq.crick

    # Watson 3' overhang
    watson_ovhg = seq.watson_ovhg()
    if watson_ovhg < 0:
        new_watson = new_watson + reverse_complement(seq.crick[:-watson_ovhg])

    # Crick 3' overhang
    elif watson_ovhg > 0:
        new_crick = reverse_complement(seq.watson[-watson_ovhg:]) + new_crick

    return _Dseq(new_watson, new_crick, seq.ovhg)


def fill_dseq(seq: _Dseq):
    return fill_left(fill_right(seq))


def reverse_complement_assembly(
    assembly: list[tuple[int, int, Location, Location]], fragments: list[_Dseqrecord]
) -> list[tuple[int, int, Location, Location]]:
    """Complement an assembly, i.e. reverse the order of the fragments and the orientation of the overlaps."""
    new_assembly = list()
    for u, v, locu, locv in assembly:
        f_u = fragments[abs(u) - 1]
        f_v = fragments[abs(v) - 1]
        new_assembly.append((-v, -u, locv._flip(len(f_v)), locu._flip(len(f_u))))
    return new_assembly[::-1]


def filter_linear_subassemblies(linear_assemblies, circular_assemblies, fragments):
    """Remove linear assemblies which are sub-assemblies of circular assemblies"""
    all_circular_assemblies = circular_assemblies + [
        reverse_complement_assembly(c, fragments) for c in circular_assemblies
    ]
    filtered_assemblies = [
        assem for assem in linear_assemblies if not any(is_sublist(assem, c, True) for c in all_circular_assemblies)
    ]
    # I don't think the line below is necessary, but just in case
    # filtered_assemblies = [l for l in filtered_assemblies if not any(is_sublist(reverse_complement_assembly(l, fragments), c, True) for c in all_circular_assemblies)]
    return filtered_assemblies


def remove_subassemblies(assemblies):
    """Filter out subassemblies, i.e. assemblies that are contained within another assembly.

    For example:
        [(1, 2, '1[8:14]:2[1:7]'), (2, 3, '2[10:17]:3[1:8]')]
        [(1, 2, '1[8:14]:2[1:7]')]
    The second one is a subassembly of the first one.
    """

    # Sort by length, longest first
    assemblies = sorted(assemblies, key=len, reverse=True)

    filtered_assemblies = list()
    for assembly in assemblies:
        # Check if this assembly is a subassembly of any of the assemblies we have already found
        if not any(is_sublist(assembly, a) for a in filtered_assemblies):
            filtered_assemblies.append(assembly)

    return filtered_assemblies


def assembly2str(assembly):
    """Convert an assembly to a string representation, for example:
    ((1, 2, [8:14], [1:7]),(2, 3, [10:17], [1:8]))
    becomes:
    ('1[8:14]:2[1:7]', '2[10:17]:3[1:8]')

    The reason for this is that by default, a feature '[8:14]' when present in a tuple
    is printed to the console as `SimpleLocation(ExactPosition(8), ExactPosition(14), strand=1)` (very long).
    """
    return str(tuple(f'{u}{lu}:{v}{lv}' for u, v, lu, lv in assembly))


def assembly2str_tuple(assembly):
    """Convert an assembly to a string representation, like
    ((1, 2, [8:14], [1:7]),(2, 3, [10:17], [1:8]))
    """
    return str(tuple((u, v, str(lu), str(lv)) for u, v, lu, lv in assembly))


def assembly_has_mismatches(fragments, assembly):
    for u, v, loc_u, loc_v in assembly:
        seq_u = fragments[u - 1] if u > 0 else fragments[-u - 1].reverse_complement()
        seq_v = fragments[v - 1] if v > 0 else fragments[-v - 1].reverse_complement()
        # TODO: Check issue where extraction failed, and whether it would give problems here
        if str(loc_u.extract(seq_u).seq).upper() != str(loc_v.extract(seq_v).seq).upper():
            return True
    return False


def assembly_is_circular(assembly, fragments):
    """
    Note: This does not work for insertion assemblies, that's why assemble takes the optional argument is_insertion.
    """
    if assembly[0][0] != assembly[-1][1]:
        return False
    elif isinstance(fragments[abs(assembly[0][0]) - 1], _Dseqrecord) and fragments[abs(assembly[0][0]) - 1].circular:
        return True
    else:
        return _location_boundaries(assembly[0][2])[0] > _location_boundaries(assembly[-1][3])[0]


def assemble(fragments, assembly, is_insertion=False):
    """Execute an assembly, from the representation returned by get_linear_assemblies or get_circular_assemblies."""

    if is_insertion:
        is_circular = False
    else:
        is_circular = assembly_is_circular(assembly, fragments)

    subfragment_representation = edge_representation2subfragment_representation(assembly, is_circular)

    # Sanity check
    for asm_edge in assembly:
        u, v, loc_u, loc_v = asm_edge
        f_u = fragments[u - 1] if u > 0 else fragments[-u - 1].reverse_complement()
        f_v = fragments[v - 1] if v > 0 else fragments[-v - 1].reverse_complement()
        seq_u = str(loc_u.extract(f_u).seq).upper()
        seq_v = str(loc_v.extract(f_v).seq).upper()
        if seq_u != seq_v:
            raise ValueError('Mismatch in assembly')

    # We transform into Dseqrecords (for primers)
    dseqr_fragments = [f if isinstance(f, _Dseqrecord) else _Dseqrecord(f) for f in fragments]
    subfragments = get_assembly_subfragments(dseqr_fragments, subfragment_representation)

    # Length of the overlaps between consecutive assembly fragments
    fragment_overlaps = [len(e[-1]) for e in assembly]

    out_dseqrecord = _Dseqrecord(subfragments[0])

    for fragment, overlap in zip(subfragments[1:], fragment_overlaps):
        # Shift the features of the right fragment to the left by `overlap`
        new_features = [f._shift(len(out_dseqrecord) - overlap) for f in fragment.features]
        # Join the left sequence including the overlap with the right sequence without the overlap
        # we use fill_right / fill_left so that it works for ligation of sticky ends
        out_dseqrecord = _Dseqrecord(
            fill_right(out_dseqrecord.seq) + fill_left(fragment.seq)[overlap:],
            features=out_dseqrecord.features + new_features,
        )

    # For circular assemblies, close the loop and wrap origin-spanning features
    if is_circular:
        overlap = fragment_overlaps[-1]

        # Special case for blunt circularisation
        if overlap == 0:
            return out_dseqrecord.looped()

        # Remove trailing overlap
        out_dseqrecord = _Dseqrecord(
            fill_dseq(out_dseqrecord.seq)[:-overlap], features=out_dseqrecord.features, circular=True
        )
        for feature in out_dseqrecord.features:
            start, end = _location_boundaries(feature.location)
            if start >= len(out_dseqrecord) or end > len(out_dseqrecord):
                # Wrap around the origin
                feature.location = _shift_location(feature.location, 0, len(out_dseqrecord))

    return out_dseqrecord


def annotate_primer_binding_sites(
    input_dseqr: _Dseqrecord, fragments: list[_Dseqrecord], assembly: list[tuple[int, int, Location, Location]]
) -> _Dseqrecord:
    """Annotate the primer binding sites in a Dseqrecord."""
    fwd, _, rvs = fragments
    start_rvs = len(input_dseqr) - len(rvs)

    output_dseqr = copy.deepcopy(input_dseqr)
    output_dseqr.add_feature(
        x=0, y=len(fwd), type_='primer_bind', strand=1, label=[fwd.name], note=['sequence: ' + str(fwd.seq)]
    )
    output_dseqr.add_feature(
        x=start_rvs,
        y=len(output_dseqr),
        type_='primer_bind',
        strand=-1,
        label=[rvs.name],
        note=['sequence: ' + str(rvs.seq)],
    )
    return output_dseqr


def edge_representation2subfragment_representation(assembly, is_circular):
    """
    Turn this kind of edge representation fragment 1, fragment 2, right edge on 1, left edge on 2
    a = [(1, 2, 'loc1a', 'loc2a'), (2, 3, 'loc2b', 'loc3b'), (3, 1, 'loc3c', 'loc1c')]
    Into this: fragment 1, left edge on 1, right edge on 1
    b = [(1, 'loc1c', 'loc1a'), (2, 'loc2a', 'loc2b'), (3, 'loc3b', 'loc3c')]
    """

    if is_circular:
        temp = list(assembly[-1:]) + list(assembly)
    else:
        temp = [(None, assembly[0][0], None, None)] + list(assembly) + [(assembly[-1][1], None, None, None)]
    edge_pairs = zip(temp, temp[1:])
    subfragment_representation = list()
    for (_u1, v1, _, start_location), (_u2, _v2, end_location, _) in edge_pairs:
        subfragment_representation.append((v1, start_location, end_location))

    return tuple(subfragment_representation)


def subfragment_representation2edge_representation(assembly, is_circular):
    """
    Turn this kind of subfragment representation fragment 1, left edge on 1, right edge on 1
    a = [(1, 'loc1c', 'loc1a'), (2, 'loc2a', 'loc2b'), (3, 'loc3b', 'loc3c')]
    Into this: fragment 1, fragment 2, right edge on 1, left edge on 2
    b = [(1, 2, 'loc1a', 'loc2a'), (2, 3, 'loc2b' 'loc3b'), (3, 1, 'loc3c', 'loc1c')]
    """

    edge_representation = []

    # Iterate through the assembly pairwise to create the edge representation
    for i in range(len(assembly) - 1):
        frag1, left1, right1 = assembly[i]
        frag2, left2, right2 = assembly[i + 1]
        # Create the edge between the current and next fragment
        edge_representation.append((frag1, frag2, right1, left2))

    if is_circular:
        # Add the edge from the last fragment back to the first
        frag_last, left_last, right_last = assembly[-1]
        frag_first, left_first, right_first = assembly[0]
        edge_representation.append((frag_last, frag_first, right_last, left_first))

    return tuple(edge_representation)


def get_assembly_subfragments(fragments: list[_Dseqrecord], subfragment_representation):
    """From the fragment representation returned by edge_representation2subfragment_representation, get the subfragments that are joined together.

    Subfragments are the slices of the fragments that are joined together

    For example:
    ```
      --A--
    TACGTAAT
      --B--
     TCGTAACGA

    Gives: TACGTAA / CGTAACGA
    ```
    To reproduce:
    ```
    a = Dseqrecord('TACGTAAT')
    b = Dseqrecord('TCGTAACGA')
    f = Assembly([a, b], limit=5)
    a0 = f.get_linear_assemblies()[0]
    print(assembly2str(a0))
    a0_subfragment_rep =edge_representation2subfragment_representation(a0, False)
    for f in get_assembly_subfragments([a, b], a0_subfragment_rep):
        print(f.seq)

    # prints TACGTAA and CGTAACGA
    ```

    Subfragments: `cccccgtatcgtgt`, `atcgtgtactgtcatattc`
    """
    subfragments = list()
    for node, start_location, end_location in subfragment_representation:
        seq = fragments[node - 1] if node > 0 else fragments[-node - 1].reverse_complement()
        subfragments.append(extract_subfragment(seq, start_location, end_location))
    return subfragments


def extract_subfragment(seq: _Dseqrecord, start_location: Location, end_location: Location):
    """Extract a subfragment from a sequence, given the start and end locations of the subfragment."""
    start = 0 if start_location is None else _location_boundaries(start_location)[0]
    end = None if end_location is None else _location_boundaries(end_location)[1]

    # Special case, some of it could be handled by better Dseqrecord slicing in the future
    if (
        seq.circular
        and start_location is not None
        and end_location is not None
        and _locations_overlap(start_location, end_location, len(seq))
    ):
        # The overhang is different for origin-spanning features, for instance
        # for a feature join{[12:13], [0:3]} in a sequence of length 13, the overhang
        # is -4, not 9
        ovhg = start - end if end > start else start - end - len(seq)
        # edge case
        if abs(ovhg) == len(seq):
            ovhg = 0
        dummy_cut = ((start, ovhg), None)
        open_seq = seq.apply_cut(dummy_cut, dummy_cut)
        return _Dseqrecord(fill_dseq(open_seq.seq), features=open_seq.features)

    return seq[start:end]


def is_sublist(sublist, my_list, my_list_is_cyclic=False):
    """Returns True if sublist is a sublist of my_list (can be treated as cyclic), False otherwise.

    Examples
    --------
    >>> is_sublist([1, 2], [1, 2, 3], False)
    True
    >>> is_sublist([1, 2], [1, 3, 2], False)
    False

    # See the case here for cyclic lists
    >>> is_sublist([3, 1], [1, 2, 3], False)
    False
    >>> is_sublist([3, 1], [1, 2, 3], True)
    True
    """
    n = len(sublist)
    if my_list_is_cyclic:
        my_list = my_list + my_list
    for i in range(len(my_list) - n + 1):
        # Just in case tuples were passed
        if list(my_list[i : i + n]) == list(sublist):
            return True
    return False


def circular_permutation_min_abs(lst):
    """Returns the circular permutation of lst with the smallest absolute value first.

    Examples
    --------
    >>> circular_permutation_min_abs([1, 2, 3])
    [1, 2, 3]
    >>> circular_permutation_min_abs([3, 1, 2])
    [1, 2, 3]
    """
    min_abs_index = min(range(len(lst)), key=lambda i: abs(lst[i]))
    return lst[min_abs_index:] + lst[:min_abs_index]


class Assembly:
    """Assembly of a list of linear DNA fragments into linear or circular
    constructs. The Assembly is meant to replace the Assembly method as it
    is easier to use. Accepts a list of Dseqrecords (source fragments) to
    initiate an Assembly object. Several methods are available for analysis
    of overlapping sequences, graph construction and assembly.

    The assembly contains a directed graph, where nodes represent fragments and
    edges represent overlaps between fragments. :
    - The node keys are integers, representing the index of the fragment in the
    input list of fragments. The sign of the node key represents the orientation
    of the fragment, positive for forward orientation, negative for reverse orientation.
    - The edges contain the locations of the overlaps in the fragments. For an edge (u, v, key):
        - u and v are the nodes connected by the edge.
        - key is a string that represents the location of the overlap. In the format:
        'u[start:end](strand):v[start:end](strand)'.
        - Edges have a 'locations' attribute, which is a list of two FeatureLocation objects,
        representing the location of the overlap in the u and v fragment, respectively.
        - You can think of an edge as a representation of the join of two fragments.

    If fragment 1 and 2 share a subsequence of 6bp, [8:14] in fragment 1 and [1:7] in fragment 2,
    there will be 4 edges representing that overlap in the graph, for all possible
    orientations of the fragments (see add_edges_from_match for details):
    - `(1, 2, '1[8:14]:2[1:7]')`
    - `(2, 1, '2[1:7]:1[8:14]')`
    - `(-1, -2, '-1[0:6]:-2[10:16]')`
    - `(-2, -1, '-2[10:16]:-1[0:6]')`

    An assembly can be thought of as a tuple of graph edges, but instead of representing them with node indexes and keys, we represent them
    as u, v, locu, locv, where u and v are the nodes connected by the edge, and locu and locv are the locations of the overlap in the first
    and second fragment. Assemblies are then represented as:
    - Linear: ((1, 2, [8:14], [1:7]), (2, 3, [10:17], [1:8]))
    - Circular: ((1, 2, [8:14], [1:7]), (2, 3, [10:17], [1:8]), (3, 1, [12:17], [1:6]))
    Note that the first and last fragment are the same in a circular assembly.

    The following constrains are applied to remove duplicate assemblies:
    - Circular assemblies: the first subfragment is not reversed, and has the smallest index in the input fragment list.
      use_fragment_order is ignored.
    - Linear assemblies:
        - Using uid (see add_edges_from_match) to identify unique edges.

    Parameters
    ----------

    fragments : list
        a list of Dseqrecord objects.
    limit : int, optional
        The shortest shared homology to be considered
    algorithm : function, optional
        The algorithm used to determine the shared sequences.
    use_fragment_order : bool, optional
        Legacy pydna behaviour: only assemblies that start with the first fragment and end with the last are considered.
    use_all_fragments : bool, optional
        Constrain the assembly to use all fragments.

    Examples
    --------

    from assembly2 import Assembly, example_fragments, assembly2str

    asm = Assembly(example_fragments, limit=5, use_fragment_order=False)
    print('Linear ===============')
    for assembly in asm.get_linear_assemblies():
        print(' ', assembly2str(assembly))
    print('Circular =============')
    for assembly in asm.get_circular_assemblies():
        print(' ', assembly2str(assembly))

    # Prints
    Linear ===============
        ('1[8:14]:2[1:7]', '2[10:17]:3[1:8]')
        ('2[10:17]:3[1:8]', '3[12:17]:1[1:6]')
        ('3[12:17]:1[1:6]', '1[8:14]:2[1:7]')
        ('1[1:6]:3[12:17]',)
        ('2[1:7]:1[8:14]',)
        ('3[1:8]:2[10:17]',)
    Circular =============
        ('1[8:14]:2[1:7]', '2[10:17]:3[1:8]', '3[12:17]:1[1:6]')

    """

    def __init__(
        self,
        frags: list[_Dseqrecord],
        limit=25,
        algorithm=common_sub_strings,
        use_fragment_order=True,
        use_all_fragments=False,
    ):
        # TODO: allow for the same fragment to be included more than once?
        self.G = _nx.MultiDiGraph()
        # Add positive and negative nodes for forward and reverse fragments
        self.G.add_nodes_from((i + 1, {'seq': f}) for (i, f) in enumerate(frags))
        self.G.add_nodes_from((-(i + 1), {'seq': f.reverse_complement()}) for (i, f) in enumerate(frags))

        # Iterate over all possible combinations of fragments
        fragment_pairs = _itertools.combinations(filter(lambda x: x > 0, self.G.nodes), 2)
        for i, j in fragment_pairs:
            # All the relative orientations of the fragments in the pair
            for u, v in _itertools.product([i, -i], [j, -j]):
                u_seq = self.G.nodes[u]['seq']
                v_seq = self.G.nodes[v]['seq']
                matches = algorithm(u_seq, v_seq, limit)
                for match in matches:
                    self.add_edges_from_match(match, u, v, u_seq, v_seq)

        self.fragments = frags
        self.limit = limit
        self.algorithm = algorithm
        self.use_fragment_order = use_fragment_order
        self.use_all_fragments = use_all_fragments

        return

    @classmethod
    def assembly_is_valid(
        cls, fragments: list[_Dseqrecord | _Primer], assembly, is_circular, use_all_fragments, is_insertion=False
    ):
        """Function used to filter paths returned from the graph, see conditions tested below."""
        if is_circular is None:
            return False

        # Linear assemblies may get begin-1-end, begin-2-end, these are removed here.
        if len(assembly) == 0:
            return False

        if use_all_fragments and len(fragments) != len(set(flatten(map(abs, e[:2]) for e in assembly))):
            return False

        # Here we check whether subsequent pairs of fragments are compatible, for instance:
        # Compatible (overlap of 1 and 2 occurs before overlap of 2 and 3):
        # (1,2,[2:9],[0:7]), (2,3,[12:19],[0:7])
        #    -- A --
        # 1 gtatcgtgt     -- B --
        # 2   atcgtgtactgtcatattc
        # 3               catattcaa
        # Incompatible (overlap of 1 and 2 occurs after overlap of 2 and 3):
        # (1,2,[2:9],[13:20]), (2,3,[0:7],[0:7])
        #                 -- A --
        #  1 -- B --    gtatcgtgt
        #  2 catattcccccccatcgtgtactgt
        #  3 catattcaa
        # Redundant: overlap of 1 and 2 ends at the same spot as overlap of 2 and 3
        # (1,2,[2:9],[1:8]), (2,3,[0:8],[0:8])
        #    -- A --
        #  gtatcgtgt
        #   catcgtgtactgtcatattc
        #   catcgtgtactgtcatattc
        #   -- B ---
        if is_circular:
            # In a circular assembly, first and last fragment must be the same
            if assembly[0][0] != assembly[-1][1]:
                return False
            edge_pairs = zip(assembly, assembly[1:] + assembly[:1])
        else:
            edge_pairs = zip(assembly, assembly[1:])

        for (_u1, v1, _, start_location), (_u2, _v2, end_location, _) in edge_pairs:
            # Incompatible as described in figure above
            fragment = fragments[abs(v1) - 1]
            if (isinstance(fragment, _Primer) or not fragment.circular) and _location_boundaries(start_location)[
                1
            ] >= _location_boundaries(end_location)[1]:
                return False

        # Fragments are used only once
        nodes_used = [
            f[0] for f in edge_representation2subfragment_representation(assembly, is_circular or is_insertion)
        ]
        if len(nodes_used) != len(set(map(abs, nodes_used))):
            return False

        return True

    def add_edges_from_match(self, match, u: int, v: int, first: _Dseqrecord, secnd: _Dseqrecord):
        """Add edges to the graph from a match returned by an `algorithm` function (see pydna.common_substrings). For
        format of edges (see documentation of the Assembly class).

        Matches are directional, because not all `algorithm` functions return the same match for (u,v) and (v,u). For example,
        homologous recombination does but sticky end ligation does not. The function returns two edges:
        - Fragments in the orientation they were passed, with locations of the match (u, v, loc_u, loc_v)
        - Reverse complement of the fragments with inverted order, with flipped locations (-v, -u, flip(loc_v), flip(loc_u))/

        """
        x_start, y_start, length = match
        if length == 0:
            # Edge case, blunt ligation
            locs = [SimpleLocation(x_start, x_start), SimpleLocation(y_start, y_start)]
        else:
            # We use shift_location with 0 to wrap origin-spanning features
            locs = [
                _shift_location(SimpleLocation(x_start, x_start + length), 0, len(first)),
                _shift_location(SimpleLocation(y_start, y_start + length), 0, len(secnd)),
            ]

        rc_locs = [locs[0]._flip(len(first)), locs[1]._flip(len(secnd))]

        # Unique id that identifies the edge in either orientation
        uid = f'{u}{locs[0]}:{v}{locs[1]}'

        combinations = (
            (u, v, locs),
            (-v, -u, rc_locs[::-1]),
        )

        for u, v, l in combinations:
            self.G.add_edge(u, v, f'{u}{l[0]}:{v}{l[1]}', locations=l, uid=uid)

    def format_assembly_edge(self, assembly_edge):
        """Go from the (u, v, key) to the (u, v, locu, locv) format."""
        u, v, key = assembly_edge
        locu, locv = self.G.get_edge_data(u, v, key)['locations']
        return u, v, locu, locv

    def get_linear_assemblies(self, only_adjacent_edges: bool = False, max_assemblies: int = 50):
        """Get linear assemblies, applying the constrains described in __init__, ensuring that paths represent
        real assemblies (see assembly_is_valid). Subassemblies are removed (see remove_subassemblies)."""

        # Copy the graph since we will add the begin and end mock nodes
        G = _nx.MultiDiGraph(self.G)
        G.add_nodes_from(['begin', 'end'])

        if self.use_fragment_order:
            # Path must start with the first fragment and end with the last
            G.add_edge('begin', 1)
            G.add_edge('begin', -1)
            G.add_edge(len(self.fragments), 'end')
            G.add_edge(-len(self.fragments), 'end')
        else:
            for node in filter(lambda x: type(x) is int, G.nodes):
                G.add_edge('begin', node)
                G.add_edge(node, 'end')

        unique_linear_paths = self.get_unique_linear_paths(G, max_assemblies)
        possible_assemblies = self.get_possible_assembly_number(unique_linear_paths)
        if possible_assemblies > max_assemblies:
            raise ValueError(f'Too many assemblies ({possible_assemblies} pre-validation) to assemble')

        assemblies = sum(map(lambda x: self.node_path2assembly_list(x, False), unique_linear_paths), [])

        out = [a for a in assemblies if self.assembly_is_valid(self.fragments, a, False, self.use_all_fragments)]
        if only_adjacent_edges:
            out = [a for a in out if self.assembly_uses_only_adjacent_edges(a, False)]
        return remove_subassemblies(out)

    def node_path2assembly_list(self, cycle, circular: bool):
        """Convert a node path in the format [1, 2, 3] (as returned by _nx.cycles.simple_cycles) to a list of all
          possible assemblies.

        There may be multiple assemblies for a given node path, if there are several edges connecting two nodes,
        for example two overlaps between 1 and 2, and single overlap between 2 and 3 should return 3 assemblies.
        """
        combine = list()
        pairing = zip(cycle, cycle[1:] + cycle[:1]) if circular else zip(cycle, cycle[1:])
        for u, v in pairing:
            combine.append([(u, v, key) for key in self.G[u][v]])
        return [tuple(map(self.format_assembly_edge, x)) for x in _itertools.product(*combine)]

    def get_unique_linear_paths(self, G_with_begin_end: _nx.MultiDiGraph, max_paths):
        # We remove the begin and end nodes, and get all paths without edges
        # e.g. we will get [1, 2, 3] only once, even if multiple edges connect
        # 1 and 2 or 2 and 3, by converting to DiGraph.

        # Cutoff has a different meaning of what one would expect, see https://github.com/networkx/networkx/issues/2762
        node_paths = [
            x[1:-1]
            for x in limit_iterator(
                _nx.all_simple_paths(_nx.DiGraph(G_with_begin_end), 'begin', 'end', cutoff=(len(self.fragments) + 1)),
                10000,
            )
        ]

        # Remove those that contain the same node twice
        node_paths = [x for x in node_paths if len(x) == len(set(map(abs, x)))]

        if self.use_all_fragments:
            node_paths = [x for x in node_paths if len(x) == len(self.fragments)]

        # For each path, we check if there are reverse complement duplicates
        # See: https://github.com/manulera/OpenCloning_backend/issues/160
        unique_node_paths = list()
        for p in node_paths:
            if [-x for x in p[::-1]] not in unique_node_paths:
                unique_node_paths.append(p)

        return unique_node_paths

    def get_possible_assembly_number(self, paths):
        possibilities = 0
        for path in paths:
            this_path = 1
            for u, v in zip(path, path[1:]):
                if v in self.G[u]:
                    this_path *= len(self.G[u][v])
            possibilities += this_path
        return possibilities

    def get_circular_assemblies(self, only_adjacent_edges: bool = False, max_assemblies: int = 50):
        """Get circular assemblies, applying the constrains described in __init__, ensuring that paths represent
        real assemblies (see assembly_is_valid)."""
        # The constrain of circular sequence is that the first node is the fragment with the smallest index in its initial orientation,
        # this is ensured by the circular_permutation_min_abs function + the filter below
        sorted_cycles = map(
            circular_permutation_min_abs,
            limit_iterator(_nx.cycles.simple_cycles(self.G, length_bound=len(self.fragments)), 10000),
        )
        sorted_cycles = filter(lambda x: x[0] > 0, sorted_cycles)
        # cycles.simple_cycles returns lists [1,2,3] not assemblies, see self.cycle2circular_assemblies

        # We apply constrains already here because sometimes the combinatorial explosion is too large
        if self.use_all_fragments:
            sorted_cycles = [c for c in sorted_cycles if len(c) == len(self.fragments)]

        # Remove cycles with duplicates
        sorted_cycles = [c for c in sorted_cycles if len(c) == len(set(map(abs, c)))]
        possible_assembly_number = self.get_possible_assembly_number([c + c[:1] for c in sorted_cycles])
        if possible_assembly_number > max_assemblies:
            raise ValueError(f'Too many assemblies ({possible_assembly_number} pre-validation) to assemble')

        assemblies = sum(map(lambda x: self.node_path2assembly_list(x, True), sorted_cycles), [])

        out = [a for a in assemblies if self.assembly_is_valid(self.fragments, a, True, self.use_all_fragments)]
        if only_adjacent_edges:
            out = [a for a in out if self.assembly_uses_only_adjacent_edges(a, True)]
        return out

    def format_insertion_assembly(self, assembly):
        """Sorts the fragment representing a cycle so that they represent an insertion assembly if possible,
        else returns None.

        Here we check if one of the joins between fragments represents the edges of an insertion assembly
        The fragment must be linear, and the join must be as indicated below

        ```
        --------         -------           Fragment 1
            ||            ||
            xxxxxxxx      ||               Fragment 2
                  ||      ||
                  oooooooooo               Fragment 3
        ```
        The above example will be [(1, 2, [4:6], [0:2]), (2, 3, [6:8], [0:2]), (3, 1, [8:10], [9:11)])]

        These could be returned in any order by simple_cycles, so we sort the edges so that the first
        and last `u` and `v` match the fragment that gets the insertion (1 in the example above).
        """
        edge_pair_index = list()

        # Pair edges with one another
        for i, ((_u1, v1, _, start_location), (_u2, _v2, end_location, _)) in enumerate(
            zip(assembly, assembly[1:] + assembly[:1])
        ):
            fragment = self.fragments[abs(v1) - 1]
            # Find the pair of edges that should be last and first  ((3, 1, [8:10], [9:11)]), (1, 2, [4:6], [0:2]) in
            # the example above. Only one of the pairs of edges should satisfy this condition for the topology to make sense.
            right_of_insertion = _location_boundaries(start_location)[0]
            left_of_insertion = _location_boundaries(end_location)[0]
            if not fragment.circular and (
                right_of_insertion > left_of_insertion
                or _locations_overlap(start_location, end_location, len(fragment))
            ):
                edge_pair_index.append(i)

        if len(edge_pair_index) != 1:
            return None

        shift_by = (edge_pair_index[0] + 1) % len(assembly)
        return assembly[shift_by:] + assembly[:shift_by]

    def get_insertion_assemblies(self, only_adjacent_edges: bool = False, max_assemblies: int = 50):
        """Assemblies that represent the insertion of a fragment or series of fragment inside a linear construct. For instance,
        digesting CCCCGAATTCCCCGAATTC with EcoRI and inserting the fragment with two overhangs into the EcoRI site of AAAGAATTCAAA.
        This is not so much meant for the use-case of linear fragments that represent actual linear fragments, but for linear
        fragments that represent a genome region. This can then be used to simulate homologous recombination.
        """
        if only_adjacent_edges:
            raise NotImplementedError('only_adjacent_edges not implemented for insertion assemblies')

        cycles = limit_iterator(_nx.cycles.simple_cycles(self.G), 10000)

        # We apply constrains already here because sometimes the combinatorial explosion is too large
        if self.use_all_fragments:
            cycles = [c for c in cycles if len(c) == len(self.fragments)]

        # Remove cycles with duplicates
        cycles = [c for c in cycles if len(c) == len(set(map(abs, c)))]

        possible_assembly_number = self.get_possible_assembly_number([c + c[:1] for c in cycles])

        if possible_assembly_number > max_assemblies:
            raise ValueError(f'Too many assemblies ({possible_assembly_number} pre-validation) to assemble')

        # We find cycles first
        iterator = limit_iterator(_nx.cycles.simple_cycles(self.G), 10000)
        assemblies = sum(map(lambda x: self.node_path2assembly_list(x, True), iterator), [])
        # We select those that contain exactly only one suitable edge
        assemblies = [b for a in assemblies if (b := self.format_insertion_assembly(a)) is not None]
        # First fragment should be in the + orientation
        assemblies = list(filter(lambda x: x[0][0] > 0, assemblies))
        return [
            a
            for a in assemblies
            if self.assembly_is_valid(self.fragments, a, False, self.use_all_fragments, is_insertion=True)
        ]

    def assemble_linear(self, only_adjacent_edges: bool = False, max_assemblies: int = 50):
        """Assemble linear constructs, from assemblies returned by self.get_linear_assemblies."""
        assemblies = self.get_linear_assemblies(only_adjacent_edges, max_assemblies)
        return [assemble(self.fragments, a) for a in assemblies]

    def assemble_circular(self, only_adjacent_edges: bool = False, max_assemblies: int = 50):
        """Assemble circular constructs, from assemblies returned by self.get_circular_assemblies."""
        assemblies = self.get_circular_assemblies(only_adjacent_edges, max_assemblies)
        return [assemble(self.fragments, a) for a in assemblies]

    def assemble_insertion(self, only_adjacent_edges: bool = False):
        """Assemble insertion constructs, from assemblies returned by self.get_insertion_assemblies."""
        assemblies = self.get_insertion_assemblies(only_adjacent_edges)
        return [assemble(self.fragments, a, is_insertion=True) for a in assemblies]

    def get_locations_on_fragments(self) -> dict[int, dict[str, list[Location]]]:
        """Get a dictionary where the keys are the nodes in the graph, and the values are dictionaries with keys
        `left`, `right`, containing (for each fragment) the locations where the fragment is joined to another fragment on its left
        and right side. The values in `left` and `right` are often the same, except in restriction-ligation with partial overlap enabled,
        where we can end up with a situation like this:

        GGTCTCCCCAATT and aGGTCTCCAACCAA as fragments

        # Partial overlap in assembly 1[9:11]:2[8:10]
        GGTCTCCxxAACCAA
        CCAGAGGGGTTxxTT

        # Partial overlap in 2[10:12]:1[7:9]
        aGGTCTCCxxCCAATT
        tCCAGAGGTTGGxxAA

        Would return
        {
            1: {'left': [7:9], 'right': [9:11]},
            2: {'left': [8:10], 'right': [10:12]},
            -1: {'left': [2:4], 'right': [4:6]},
            -2: {'left': [2:4], 'right': [4:6]}
        }

        """

        locations_on_fragments = dict()
        for node in self.G.nodes:
            this_dict = {'left': list(), 'right': list()}
            for edge in self.G.edges(data=True):
                for i, key in enumerate(['right', 'left']):
                    if edge[i] == node:
                        edge_location = edge[2]['locations'][i]
                        if edge_location not in this_dict[key]:
                            this_dict[key].append(edge_location)
            this_dict['left'] = sorted(this_dict['left'], key=lambda x: _location_boundaries(x)[0])
            this_dict['right'] = sorted(this_dict['right'], key=lambda x: _location_boundaries(x)[0])
            locations_on_fragments[node] = this_dict

        return locations_on_fragments

    def assembly_uses_only_adjacent_edges(self, assembly, is_circular: bool) -> bool:
        """
        Check whether only adjacent edges within each fragment are used in the assembly. This is useful to check if a cut and ligate assembly is valid,
          and prevent including partially digested fragments. For example, imagine the following fragment being an input for a digestion
          and ligation assembly, where the enzyme cuts at the sites indicated by the vertical lines:

                 x       y       z
          -------|-------|-------|---------

          We would only want assemblies that contain subfragments start-x, x-y, y-z, z-end, and not start-x, y-end, for instance.
          The latter would indicate that the fragment was partially digested.
        """

        locations_on_fragments = self.get_locations_on_fragments()
        for node in locations_on_fragments:
            fragment_len = len(self.fragments[abs(node) - 1])
            for side in ['left', 'right']:
                locations_on_fragments[node][side] = gather_overlapping_locations(
                    locations_on_fragments[node][side], fragment_len
                )

        allowed_location_pairs = dict()
        for node in locations_on_fragments:
            if not is_circular:
                # We add the existing ends of the fragment
                left = [(None,)] + locations_on_fragments[node]['left']
                right = locations_on_fragments[node]['right'] + [(None,)]

            else:
                # For circular assemblies, we add the first location at the end
                # to allow for the last edge to be used
                left = locations_on_fragments[node]['left']
                right = locations_on_fragments[node]['right'][1:] + locations_on_fragments[node]['right'][:1]

            pairs = list()
            for pair in zip(left, right):
                pairs += list(_itertools.product(*pair))
            allowed_location_pairs[node] = pairs

        fragment_assembly = edge_representation2subfragment_representation(assembly, is_circular)
        for node, start_location, end_location in fragment_assembly:
            if (start_location, end_location) not in allowed_location_pairs[node]:
                return False
        return True

    def __repr__(self):
        # https://pyformat.info
        return _pretty_str(
            'Assembly\n'
            'fragments..: {sequences}\n'
            'limit(bp)..: {limit}\n'
            'G.nodes....: {nodes}\n'
            'algorithm..: {al}'.format(
                sequences=' '.join('{}bp'.format(len(x)) for x in self.fragments),
                limit=self.limit,
                nodes=self.G.order(),
                al=self.algorithm.__name__,
            )
        )


class PCRAssembly(Assembly):
    def __init__(self, frags: list[_Dseqrecord | _Primer], limit=25, mismatches=0):

        value_error = ValueError(
            'PCRAssembly assembly must be initialised with a list/tuple of primer, template, primer'
        )
        if len(frags) != 3:
            raise value_error

        # Validate the inputs: should be a series of primer, template, primer
        wrong_fragment_class = (
            not isinstance(frags[0], _Primer),
            isinstance(frags[1], _Primer),
            not isinstance(frags[2], _Primer),
        )
        if any(wrong_fragment_class):
            raise value_error

        # TODO: allow for the same fragment to be included more than once?
        self.G = _nx.MultiDiGraph()
        # Add positive and negative nodes for forward and reverse fragments
        self.G.add_nodes_from((i + 1, {'seq': f}) for (i, f) in enumerate(frags))
        self.G.add_nodes_from((-(i + 1), {'seq': f.reverse_complement()}) for (i, f) in enumerate(frags))

        pairs = list()
        primer_ids = list()
        for i in range(0, len(frags), 3):
            # primer, template, primer
            p1, t, p2 = (i + 1, i + 2, i + 3)
            primer_ids += [p1, p2]
            pairs += list(_itertools.product([p1, p2], [t, -t]))
            pairs += list(_itertools.product([t, -t], [-p1, -p2]))

        for u, v in pairs:
            u_seq = self.G.nodes[u]['seq']
            v_seq = self.G.nodes[v]['seq']
            matches = alignment_sub_strings(u_seq, v_seq, limit, mismatches)
            for match in matches:
                self.add_edges_from_match(match, u, v, u_seq, v_seq)

        # These two are constrained
        self.use_fragment_order = False
        self.use_all_fragments = True

        self.fragments = frags
        self.limit = limit
        self.algorithm = alignment_sub_strings

        return

    def get_linear_assemblies(self, only_adjacent_edges: bool = False, max_assemblies: int = 50):
        if only_adjacent_edges:
            raise NotImplementedError('only_adjacent_edges not implemented for PCR assemblies')

        return super().get_linear_assemblies(max_assemblies=max_assemblies)

    def get_circular_assemblies(self, only_adjacent_edges: bool = False):
        raise NotImplementedError('get_circular_assemblies not implemented for PCR assemblies')

    def get_insertion_assemblies(self, only_adjacent_edges: bool = False):
        raise NotImplementedError('get_insertion_assemblies not implemented for PCR assemblies')


class SingleFragmentAssembly(Assembly):
    """
    An assembly that represents the circularisation or splicing of a single fragment.
    """

    def __init__(self, frags: [_Dseqrecord], limit=25, algorithm=common_sub_strings):

        if len(frags) != 1:
            raise ValueError('SingleFragmentAssembly assembly must be initialised with a single fragment')
        # TODO: allow for the same fragment to be included more than once?
        self.G = _nx.MultiDiGraph()
        frag = frags[0]
        # Add positive and negative nodes for forward and reverse fragments
        self.G.add_node(1, seq=frag)

        matches = algorithm(frag, frag, limit)
        for match in matches:
            self.add_edges_from_match(match, 1, 1, frag, frag)

        # To avoid duplicated outputs
        self.G.remove_edges_from([(-1, -1)])

        # These two are constrained
        self.use_fragment_order = True
        self.use_all_fragments = True

        self.fragments = frags
        self.limit = limit
        self.algorithm = algorithm

        return

    def get_circular_assemblies(self, only_adjacent_edges: bool = False, max_assemblies: int = 50):
        # We don't want the same location twice
        assemblies = filter(
            lambda x: x[0][2] != x[0][3], super().get_circular_assemblies(only_adjacent_edges, max_assemblies)
        )
        return [a for a in assemblies if self.assembly_is_valid(self.fragments, a, True, self.use_all_fragments)]

    def get_insertion_assemblies(self, only_adjacent_edges: bool = False, max_assemblies: int = 50):
        """This could be renamed splicing assembly, but the essence is similar"""

        if only_adjacent_edges:
            raise NotImplementedError('only_adjacent_edges not implemented for insertion assemblies')

        def splicing_assembly_filter(x):
            # We don't want the same location twice
            if x[0][2] == x[0][3]:
                return False
            # We don't want to get overlap only (e.g. GAATTCcatGAATTC giving GAATTC)
            left_start, _ = _location_boundaries(x[0][2])
            _, right_end = _location_boundaries(x[0][3])
            if left_start == 0 and right_end == len(self.fragments[0]):
                return False
            return True

        # We don't want the same location twice
        assemblies = filter(splicing_assembly_filter, super().get_insertion_assemblies(max_assemblies=max_assemblies))
        return [
            a
            for a in assemblies
            if self.assembly_is_valid(self.fragments, a, False, self.use_all_fragments, is_insertion=True)
        ]

    def get_linear_assemblies(self):
        raise NotImplementedError('Linear assembly does not make sense')


example_fragments = (
    _Dseqrecord('AacgatCAtgctcc', name='a'),
    _Dseqrecord('TtgctccTAAattctgc', name='b'),
    _Dseqrecord('CattctgcGAGGacgatG', name='c'),
)


'CattctgcGAGGacgatCAtgctcc'

linear_results = (
    _Dseqrecord('AacgatCAtgctccTAAattctgcGAGGacgatG', name='abc'),
    _Dseqrecord('ggagcaTGatcgtCCTCgcagaatG', name='ac_rc'),
    _Dseqrecord('AacgatG', name='ac'),
)


circular_results = (_Dseqrecord('acgatCAtgctccTAAattctgcGAGG', name='abc', circular=True),)
