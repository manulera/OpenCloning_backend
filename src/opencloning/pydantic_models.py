from pydantic import BaseModel, Field, model_validator, field_validator, Discriminator, Tag
from typing import Optional, List, Union, Annotated

from ._version import __version__

from Bio.SeqFeature import (
    Location,
)
from Bio.SeqRecord import SeqRecord as _SeqRecord
from opencloning_linkml.datamodel import (
    OligoHybridizationSource as _OligoHybridizationSource,
    PolymeraseExtensionSource as _PolymeraseExtensionSource,
    GenomeCoordinatesSource as _GenomeCoordinatesSource,
    RepositoryIdSource as _RepositoryIdSource,
    ManuallyTypedSource as _ManuallyTypedSource,
    UploadedFileSource as _UploadedFileSource,
    SequenceFileFormat as _SequenceFileFormat,
    RestrictionEnzymeDigestionSource as _RestrictionEnzymeDigestionSource,
    TextFileSequence as _TextFileSequence,
    AssemblySource as _AssemblySource,
    PCRSource as _PCRSource,
    HomologousRecombinationSource as _HomologousRecombinationSource,
    GibsonAssemblySource as _GibsonAssemblySource,
    RestrictionAndLigationSource as _RestrictionAndLigationSource,
    LigationSource as _LigationSource,
    CRISPRSource as _CRISPRSource,
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
    SourceInput as _SourceInput,
    OpenDNACollectionsSource as _OpenDNACollectionsSource,
    Primer as PrimerModel,
)
from pydna.assembly2 import (
    edge_representation2subfragment_representation,
    subfragment_representation2edge_representation,
)
from pydna.opencloning_models import SequenceLocationStr


SequenceFileFormat = _SequenceFileFormat


class SourceInput(_SourceInput):
    pass


# Sources =========================================


def input_discriminator(v) -> str | None:
    """
    Discriminator that yields SourceInput by default
    """
    if isinstance(v, dict):
        input_type = v.get('type', None)
        if input_type is None:
            return 'SourceInput'
        else:
            return input_type
    elif isinstance(v, SourceInput):
        return v.type
    return None


class SourceCommonClass(BaseModel):
    input: Optional[List[SourceInput]] = Field(
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
    pass


class OpenDNACollectionsSource(SourceCommonClass, _OpenDNACollectionsSource):
    pass


class SEVASource(SourceCommonClass, _SEVASource):
    pass


class GenomeCoordinatesSource(SourceCommonClass, _GenomeCoordinatesSource):
    pass


class AnnotationSource(SourceCommonClass, _AnnotationSource):
    pass


class ReverseComplementSource(SourceCommonClass, _ReverseComplementSource):
    pass


class RestrictionEnzymeDigestionSource(SourceCommonClass, _RestrictionEnzymeDigestionSource):
    pass


class AssemblyFragment(_AssemblyFragment, SourceInput):
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
    input: Optional[
        List[
            Annotated[
                Union[
                    Annotated[SourceInput, Tag('SourceInput')],
                    Annotated['AssemblyFragment', Tag('AssemblyFragment')],
                ],
                Discriminator(input_discriminator),
            ]
        ]
    ] = Field(
        default_factory=list,
        description="""The inputs to this source. If the source represents external import of a sequence, it's empty.""",
        json_schema_extra={'linkml_meta': {'alias': 'input', 'domain_of': ['Source'], 'slot_uri': 'schema:object'}},
    )

    def get_assembly_plan(self, fragments: list[_SeqRecord]) -> tuple:
        """Returns the assembly plan"""
        subf = [f.to_fragment_tuple(fragments) for f in self.input if f.type == 'AssemblyFragment']
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
            input=assembly_fragments,
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
        source = super().from_assembly(assembly, id, False, fragments)
        source.input += [SourceInput(sequence=guide) for guide in guides]
        return source


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

    def add_primer(self, primer: PrimerModel):
        if primer in self.primers:
            return
        primer.id = self.next_id()
        self.primers.append(primer)

    def next_id(self):
        return max([s.id for s in self.sources + self.sequences + self.primers], default=0) + 1

    def add_source_and_sequence(self, source: SourceCommonClass, sequence: _TextFileSequence):
        if source in self.sources:
            if sequence not in self.sequences:
                raise ValueError(
                    f"Source {source.id} already exists in the cloning strategy, but sequence {sequence.id} it's not its output."
                )
            return
        new_id = self.next_id()
        source.id = new_id
        self.sources.append(source)
        sequence.id = new_id
        self.sequences.append(sequence)

    def all_children_source_ids(self, source_id: int, source_children: list | None = None) -> list[int]:
        """Returns the ids of all source children ids of a source"""
        source = next(s for s in self.sources if s.id == source_id)
        if source_children is None:
            source_children = []

        sources_that_take_output_as_input = [s for s in self.sources if source.id in [inp.sequence for inp in s.input]]
        new_source_ids = [s.id for s in sources_that_take_output_as_input]

        source_children.extend(new_source_ids)
        for new_source_id in new_source_ids:
            self.all_children_source_ids(new_source_id, source_children)
        return source_children


class PrimerDesignQuery(BaseModel):
    model_config = {'arbitrary_types_allowed': True}
    sequence: _TextFileSequence
    location: SequenceLocationStr
    forward_orientation: bool = True

    @field_validator('location', mode='before')
    def parse_location(cls, v):
        return SequenceLocationStr.field_validator(v)
