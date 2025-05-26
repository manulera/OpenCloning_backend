from pydantic import BaseModel, Field, model_validator, field_validator
from typing import Optional, List
from pydantic_core import core_schema
from ._version import __version__

from Bio.SeqFeature import (
    SeqFeature,
    Location,
    SimpleLocation,
    FeatureLocation as BioFeatureLocation,
    LocationParserError,
)
from Bio.SeqIO.InsdcIO import _insdc_location_string as format_feature_location
from Bio.Restriction.Restriction import RestrictionType, RestrictionBatch
from Bio.SeqRecord import SeqRecord as _SeqRecord
from pydna.primer import Primer as _PydnaPrimer
from opencloning_linkml.datamodel import (
    OligoHybridizationSource as _OligoHybridizationSource,
    PolymeraseExtensionSource as _PolymeraseExtensionSource,
    GenomeCoordinatesSource as _GenomeCoordinatesSource,
    RepositoryIdSource as _RepositoryIdSource,
    ManuallyTypedSource as _ManuallyTypedSource,
    UploadedFileSource as _UploadedFileSource,
    SequenceFileFormat as _SequenceFileFormat,
    RestrictionEnzymeDigestionSource as _RestrictionEnzymeDigestionSource,
    RestrictionSequenceCut as _RestrictionSequenceCut,
    TextFileSequence as _TextFileSequence,
    AssemblySource as _AssemblySource,
    PCRSource as _PCRSource,
    HomologousRecombinationSource as _HomologousRecombinationSource,
    GibsonAssemblySource as _GibsonAssemblySource,
    RestrictionAndLigationSource as _RestrictionAndLigationSource,
    LigationSource as _LigationSource,
    CRISPRSource as _CRISPRSource,
    Primer as _Primer,
    AssemblyFragment as _AssemblyFragment,
    AddgeneIdSource as _AddgeneIdSource,
    WekWikGeneIdSource as _WekWikGeneIdSource,
    BenchlingUrlSource as _BenchlingUrlSource,
    CloningStrategy as _CloningStrategy,
    OverlapExtensionPCRLigationSource as _OverlapExtensionPCRLigationSource,
    SnapGenePlasmidSource as _SnapGenePlasmidSource,
    EuroscarfSource as _EuroscarfSource,
    GatewaySource as _GatewaySource,
    InFusionSource as _InFusionSource,
    AnnotationSource as _AnnotationSource,
    IGEMSource as _IGEMSource,
    ReverseComplementSource as _ReverseComplementSource,
    SEVASource as _SEVASource,
    CreLoxRecombinationSource as _CreLoxRecombinationSource,
    InVivoAssemblySource as _InVivoAssemblySource,
)
from .assembly2 import (
    edge_representation2subfragment_representation,
    subfragment_representation2edge_representation,
)
from pydna.utils import location_boundaries, shift_location


SequenceFileFormat = _SequenceFileFormat


class TextFileSequence(_TextFileSequence):
    pass


class PrimerModel(_Primer):
    """Called PrimerModel not to be confused with the class from pydna."""

    def to_pydna_primer(self) -> _PydnaPrimer:
        """
        Convert the PrimerModel to a pydna Primer object.

        Returns:
            _PydnaPrimer: A pydna Primer object with the same sequence and name as the PrimerModel.
        """
        return _PydnaPrimer(self.sequence, name=self.name, id=str(self.id))


class SeqFeatureModel(BaseModel):
    type: str
    qualifiers: dict[str, list[str]] = {}
    location: str

    def convert_to_seq_feature(self) -> SeqFeature:
        return SeqFeature(location=Location.fromstring(self.location), type=self.type, qualifiers=self.qualifiers)

    def read_from_seq_feature(sf: SeqFeature) -> 'SeqFeatureModel':
        return SeqFeatureModel(
            type=sf.type, qualifiers=sf.qualifiers, location=format_feature_location(sf.location, None)
        )


# Sources =========================================


class SourceCommonClass:
    input: Optional[List[int]] = Field(
        default_factory=list,
        description="""The sequences that are an input to this source. If the source represents external import of a sequence, it's empty.""",
        json_schema_extra={'linkml_meta': {'alias': 'input', 'domain_of': ['Source']}},
    )


class ManuallyTypedSource(SourceCommonClass, _ManuallyTypedSource):
    """Describes a sequence that is typed manually by the user"""

    @model_validator(mode='after')
    def validate_circularity(self):
        # Do the validation instead of printing
        if self.circular:
            assert self.overhang_crick_3prime == 0, 'Circular sequences cannot have overhangs.'
            assert self.overhang_watson_3prime == 0, 'Circular sequences cannot have overhangs.'
        return self


class UploadedFileSource(SourceCommonClass, _UploadedFileSource):
    coordinates: Optional['SequenceLocationStr'] = Field(
        default=None,
        description="""If provided, coordinates within the sequence of the file to extract a subsequence""",
        json_schema_extra={'linkml_meta': {'alias': 'coordinates', 'domain_of': ['UploadedFileSource']}},
    )

    @field_validator('coordinates', mode='before')
    def parse_coordinates(cls, v):
        if v is None:
            return None
        return SequenceLocationStr.field_validator(v)


class RepositoryIdSource(SourceCommonClass, _RepositoryIdSource):
    pass


class AddgeneIdSource(SourceCommonClass, _AddgeneIdSource):
    # TODO: add this to LinkML
    # repository_name: RepositoryName = RepositoryName('addgene')
    pass


class WekWikGeneIdSource(SourceCommonClass, _WekWikGeneIdSource):
    pass


class BenchlingUrlSource(SourceCommonClass, _BenchlingUrlSource):
    pass


class SnapGenePlasmidSource(SourceCommonClass, _SnapGenePlasmidSource):
    pass


class EuroscarfSource(SourceCommonClass, _EuroscarfSource):
    pass


class IGEMSource(SourceCommonClass, _IGEMSource):

    @model_validator(mode='after')
    def validate_repository_id(self):
        file_name = self.sequence_file_url.split('/')[-1]
        assert file_name.endswith('.gb'), 'The sequence file must be a GenBank file'
        return self


class SEVASource(SourceCommonClass, _SEVASource):
    pass


class GenomeCoordinatesSource(SourceCommonClass, _GenomeCoordinatesSource):
    pass


class AnnotationSource(SourceCommonClass, _AnnotationSource):
    pass


class ReverseComplementSource(SourceCommonClass, _ReverseComplementSource):
    pass


class RestrictionSequenceCut(_RestrictionSequenceCut):

    @classmethod
    def from_cutsite_tuple(cls, cutsite_tuple: tuple[tuple[int, int], RestrictionType]):
        cut_watson, ovhg = cutsite_tuple[0]
        enzyme = str(cutsite_tuple[1])

        return cls(
            cut_watson=cut_watson,
            overhang=ovhg,
            restriction_enzyme=enzyme,
        )

    def to_cutsite_tuple(self) -> tuple[tuple[int, int], RestrictionType]:
        restriction_enzyme = RestrictionBatch(first=[self.restriction_enzyme]).pop()
        return ((self.cut_watson, self.overhang), restriction_enzyme)


class RestrictionEnzymeDigestionSource(SourceCommonClass, _RestrictionEnzymeDigestionSource):
    """Documents a restriction enzyme digestion, and the selection of one of the fragments."""

    # TODO: maybe a better way? They have to be redefined here because
    # we have overriden the original class

    left_edge: Optional[RestrictionSequenceCut] = Field(None)
    right_edge: Optional[RestrictionSequenceCut] = Field(None)

    @classmethod
    def from_cutsites(
        cls,
        left: tuple[tuple[int, int], RestrictionType],
        right: tuple[tuple[int, int], RestrictionType],
        input: list[int],
        id: int,
    ):
        return cls(
            id=id,
            left_edge=None if left is None else RestrictionSequenceCut.from_cutsite_tuple(left),
            right_edge=None if right is None else RestrictionSequenceCut.from_cutsite_tuple(right),
            input=input,
        )

    # TODO could be made into a computed field?
    def get_enzymes(self) -> list[str]:
        """Returns the enzymes used in the digestion"""
        out = list()
        if self.left_edge is not None:
            out.append(self.left_edge.restriction_enzyme)
        if self.right_edge is not None:
            out.append(self.right_edge.restriction_enzyme)
        # Unique values, sorted the same way
        return sorted(list(set(out)), key=out.index)


class SequenceLocationStr(str):
    """A string representation of a sequence location, genbank-like."""

    # TODO: this should handle origin-spanning simple locations (splitted)
    @classmethod
    def from_biopython_location(cls, location: Location):
        return cls(format_feature_location(location, None))

    @classmethod
    def from_start_and_end(cls, start: int, end: int, seq_len: int | None = None, strand: int | None = 1):
        if end >= start:
            return cls.from_biopython_location(SimpleLocation(start, end, strand=strand))
        else:
            if seq_len is None:
                raise ValueError('Sequence length is required to handle origin-spanning simple locations')
            unwrapped_location = SimpleLocation(start, end + seq_len, strand=strand)
            wrapped_location = shift_location(unwrapped_location, 0, seq_len)
            return cls.from_biopython_location(wrapped_location)

    def to_biopython_location(self) -> BioFeatureLocation:
        return Location.fromstring(self)

    @classmethod
    def field_validator(cls, v):
        if isinstance(v, str):
            value = cls(v)
            try:
                value.to_biopython_location()
            except LocationParserError:
                raise ValueError(f'Location "{v}" is not a valid location')
            return value
        raise ValueError(f'Location must be a string or a {cls.__name__}')

    @property
    def start(self) -> int:
        return location_boundaries(self.to_biopython_location())[0]

    @property
    def end(self) -> int:
        return location_boundaries(self.to_biopython_location())[1]

    @classmethod
    def __get_pydantic_core_schema__(
        cls,
        source_type,
        handler,
    ) -> core_schema.CoreSchema:
        """Generate Pydantic core schema for SequenceLocationStr."""
        return core_schema.with_info_after_validator_function(
            cls._validate,
            core_schema.str_schema(),
        )

    @classmethod
    def _validate(cls, value: str, info):
        """Validate and create SequenceLocationStr instance."""
        return cls.field_validator(value)


class AssemblyFragment(_AssemblyFragment):
    left_location: Optional[SequenceLocationStr] = None
    right_location: Optional[SequenceLocationStr] = None

    def to_fragment_tuple(self, fragments) -> tuple[int, Location, Location]:
        fragment_ids = [int(f.id) for f in fragments]
        # By convention, these have no strand
        left_loc = None if self.left_location is None else self.left_location.to_biopython_location()
        right_loc = None if self.right_location is None else self.right_location.to_biopython_location()
        if left_loc is not None:
            left_loc.strand = None
        if right_loc is not None:
            right_loc.strand = None

        return (
            (fragment_ids.index(self.sequence) + 1) * (-1 if self.reverse_complemented else 1),
            left_loc,
            right_loc,
        )

    @field_validator('left_location', 'right_location', mode='before')
    def parse_location(cls, v):
        if v is None:
            return None
        return SequenceLocationStr.field_validator(v)


class AssemblySourceCommonClass(SourceCommonClass):
    # TODO: This is different in the LinkML model, because there it is not required,
    # and here we make it default to list.
    assembly: List[AssemblyFragment] = Field(
        default_factory=list, description="""The joins between the fragments in the assembly"""
    )

    def minimal_overlap(self):
        """Returns the minimal overlap between the fragments in the assembly"""
        all_overlaps = list()
        for f in self.assembly:
            if f.left_location is not None:
                all_overlaps.append(f.left_location.end - f.left_location.start)
            if f.right_location is not None:
                all_overlaps.append(f.right_location.end - f.right_location.start)
        return min(all_overlaps)

    def get_assembly_plan(self, fragments: list[_SeqRecord]) -> tuple:
        """Returns the assembly plan"""
        subf = [f.to_fragment_tuple(fragments) for f in self.assembly]
        return subfragment_representation2edge_representation(subf, self.circular)

    @classmethod
    def from_assembly(
        cls,
        assembly: list[tuple[int, int, Location, Location]],
        id: int,
        circular: bool,
        fragments: list[_SeqRecord],
        **kwargs,
    ):

        # Replace the positions with the actual ids
        fragment_ids = [int(f.id) for f in fragments]
        input_ids = [int(f.id) for f in fragments if not isinstance(f, _PydnaPrimer)]

        # Here the ids are still the positions in the fragments list
        fragment_assembly_positions = edge_representation2subfragment_representation(assembly, circular)
        assembly_fragments = [
            AssemblyFragment(
                sequence=fragment_ids[abs(pos) - 1],
                left_location=None if left_loc is None else SequenceLocationStr.from_biopython_location(left_loc),
                right_location=None if right_loc is None else SequenceLocationStr.from_biopython_location(right_loc),
                reverse_complemented=pos < 0,
            )
            for pos, left_loc, right_loc in fragment_assembly_positions
        ]
        return cls(
            id=id,
            input=input_ids,
            assembly=assembly_fragments,
            circular=circular,
            **kwargs,
        )


class AssemblySource(AssemblySourceCommonClass, _AssemblySource):
    pass


class PCRSource(AssemblySourceCommonClass, _PCRSource):
    pass


class LigationSource(AssemblySourceCommonClass, _LigationSource):
    pass


class HomologousRecombinationSource(AssemblySourceCommonClass, _HomologousRecombinationSource):

    # TODO: add this to LinkML
    # This can only take two inputs, the first one is the template, the second one is the insert
    # input: conlist(int, min_length=2, max_length=2)
    pass


class GibsonAssemblySource(AssemblySourceCommonClass, _GibsonAssemblySource):

    # TODO: add this to LinkML
    # input: conlist(int, min_length=1)
    pass


class OverlapExtensionPCRLigationSource(AssemblySourceCommonClass, _OverlapExtensionPCRLigationSource):
    pass


class InFusionSource(AssemblySourceCommonClass, _InFusionSource):
    pass


class InVivoAssemblySource(AssemblySourceCommonClass, _InVivoAssemblySource):
    pass


class CRISPRSource(AssemblySourceCommonClass, _CRISPRSource):

    # TODO
    # input: conlist(int, min_length=2, max_length=2)
    # circular: bool = False

    @classmethod
    def from_assembly(
        cls,
        assembly: list[tuple[int, int, Location, Location]],
        id: int,
        fragments: list[_SeqRecord],
        guides: list[int],
    ):
        return super().from_assembly(assembly, id, False, fragments, guides=guides)


class RestrictionAndLigationSource(AssemblySourceCommonClass, _RestrictionAndLigationSource):
    # TODO: add this to LinkML
    # input: conlist(int, min_length=1)

    @classmethod
    def from_assembly(
        cls,
        assembly: list[tuple[int, int, Location, Location]],
        circular: bool,
        id: int,
        fragments: list[_SeqRecord],
        restriction_enzymes=list['str'],
    ):
        return super().from_assembly(assembly, id, circular, fragments, restriction_enzymes=restriction_enzymes)


class GatewaySource(AssemblySourceCommonClass, _GatewaySource):
    @classmethod
    def from_assembly(
        cls,
        assembly: list[tuple[int, int, Location, Location]],
        circular: bool,
        id: int,
        fragments: list[_SeqRecord],
        reaction_type: str,
    ):
        return super().from_assembly(assembly, id, circular, fragments, reaction_type=reaction_type)


class CreLoxRecombinationSource(AssemblySourceCommonClass, _CreLoxRecombinationSource):
    pass


class OligoHybridizationSource(SourceCommonClass, _OligoHybridizationSource):
    pass


class PolymeraseExtensionSource(SourceCommonClass, _PolymeraseExtensionSource):
    pass


class BaseCloningStrategy(_CloningStrategy):
    # For now, we don't add anything, but the classes will not have the new methods if this is used
    # It will be used for validation for now
    primers: Optional[List[PrimerModel]] = Field(
        default_factory=list,
        description="""The primers that are used in the cloning strategy""",
        json_schema_extra={'linkml_meta': {'alias': 'primers', 'domain_of': ['CloningStrategy']}},
    )
    backend_version: Optional[str] = Field(
        default=__version__,
        description="""The version of the backend that was used to generate this cloning strategy""",
        json_schema_extra={'linkml_meta': {'alias': 'backend_version', 'domain_of': ['CloningStrategy']}},
    )

    def next_primer_id(self):
        return max([p.id for p in self.primers], default=0) + 1

    def add_primer(self, primer: PrimerModel):
        if primer in self.primers:
            return
        primer.id = self.next_primer_id()
        self.primers.append(primer)

    def next_node_id(self):
        return max([s.id for s in self.sources + self.sequences], default=0) + 1

    def add_source_and_sequence(self, source: SourceCommonClass, sequence: TextFileSequence):
        if source in self.sources:
            if sequence not in self.sequences:
                raise ValueError(
                    f"Source {source.id} already exists in the cloning strategy, but sequence {sequence.id} it's not its output."
                )
            return
        source.id = self.next_node_id()
        self.sources.append(source)
        sequence.id = self.next_node_id()
        self.sequences.append(sequence)
        source.output = sequence.id

    def all_children_source_ids(self, source_id: int, source_children: list | None = None) -> list[int]:
        """Returns the ids of all source children ids of a source"""
        source = next(s for s in self.sources if s.id == source_id)
        if source_children is None:
            source_children = []

        sources_that_take_output_as_input = [s for s in self.sources if source.output in s.input]
        new_source_ids = [s.id for s in sources_that_take_output_as_input]

        source_children.extend(new_source_ids)
        for new_source_id in new_source_ids:
            self.all_children_source_ids(new_source_id, source_children)
        return source_children


class PrimerDesignQuery(BaseModel):
    model_config = {'arbitrary_types_allowed': True}
    sequence: TextFileSequence
    location: SequenceLocationStr
    forward_orientation: bool = True

    @field_validator('location', mode='before')
    def parse_location(cls, v):
        return SequenceLocationStr.field_validator(v)
