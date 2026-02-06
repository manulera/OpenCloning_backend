import re
import networkx as nx
from Bio.Restriction.Restriction import RestrictionType
from Bio.Seq import reverse_complement
from typing import Dict, List
from pydna.dseqrecord import Dseqrecord

from pydantic import BaseModel, Field, GetJsonSchemaHandler, field_validator
from pydantic.json_schema import JsonSchemaValue
from pydantic_core import core_schema

from opencloning.endpoints.endpoint_utils import parse_restriction_enzymes

from .css_colors import CSS_COLORS


class DNASequence(str):
    """A DNA sequence string that only contains ACGT characters (uppercase)."""

    @classmethod
    def __get_pydantic_core_schema__(cls, source_type, handler) -> core_schema.CoreSchema:
        return core_schema.no_info_after_validator_function(
            cls._validate,
            core_schema.str_schema(),
            serialization=core_schema.str_schema(),
        )

    @classmethod
    def _validate(cls, value: str) -> 'DNASequence':
        if not isinstance(value, str):
            raise TypeError('DNA sequence must be a string')
        if value and not re.match(r'^[ACGT]+$', value):
            raise ValueError(f"DNA sequence must only contain uppercase ACGT characters, " f"got: {value}")
        return cls(value)

    @classmethod
    def __get_pydantic_json_schema__(cls, _core_schema, handler: GetJsonSchemaHandler) -> JsonSchemaValue:
        return {'type': 'string', 'pattern': '^[ACGT]*$'}


class Part(BaseModel):
    """Represents a single part in the syntax."""

    id: int
    name: str
    info: str = Field(default='')
    glyph: str = Field(default='')
    left_overhang: DNASequence
    right_overhang: DNASequence
    left_inside: DNASequence = Field(default='')
    right_inside: DNASequence = Field(default='')
    left_codon_start: int = Field(default=0)
    right_codon_start: int = Field(default=0)
    color: str = Field(default='')

    @property
    def key(self) -> str:
        """The key of the part is the concatenation of the left and right overhangs."""
        return f"{self.left_overhang}-{self.right_overhang}"

    @field_validator('color')
    @classmethod
    def validate_color(cls, value: str) -> str:
        """Validate color format - can be empty or valid CSS color."""
        # Color can be empty
        if not value or value.strip() == '':
            return ''

        color = value.strip()

        # Check if it's a valid CSS color name
        if color.lower() in CSS_COLORS:
            return color

        # Check hex color format (# followed by 3-6 hex digits)
        if re.match(r'^#[0-9A-Fa-f]{3,6}$', color):
            return color

        # Check rgb() format
        if re.match(r'^rgb\(\d{1,3},\d{1,3},\d{1,3}\)$', color):
            return color

        # Check rgba() format
        if re.match(r'^rgba\(\d{1,3},\d{1,3},\d{1,3},\d{1,3}\)$', color):
            return color

        # Check hsl() format
        if re.match(r'^hsl\(\d{1,3},\d{1,3},\d{1,3}\)$', color):
            return color

        # Check hsla() format
        if re.match(r'^hsla\(\d{1,3},\d{1,3},\d{1,3},\d{1,3}\)$', color):
            return color

        raise ValueError(f"Invalid color format: {color}")

    @field_validator('left_overhang', 'right_overhang')
    @classmethod
    def validate_overhang(cls, value: DNASequence) -> DNASequence:
        if len(value) != 4:
            raise ValueError(f"Overhang must be 4 characters long, got: {value}")
        return DNASequence(value)

    @field_validator('left_codon_start', 'right_codon_start')
    @classmethod
    def validate_codon_start(cls, value: int) -> int:
        if value < 0:
            raise ValueError(f"Codon start must be non-negative, got: {value}")
        return value


class Syntax(BaseModel):
    """Represents a complete syntax definition."""

    syntax_name: str = Field(alias='syntaxName')
    assembly_enzyme: str = Field(alias='assemblyEnzyme')
    domestication_enzyme: str = Field(alias='domesticationEnzyme')
    related_dois: List[str] = Field(alias='relatedDois')
    submitters: List[str]
    overhang_names: Dict[DNASequence, str] = Field(alias='overhangNames')
    parts: List[Part]

    @field_validator('submitters')
    @classmethod
    def validate_submitters(cls, value: List[str]) -> List[str]:
        for submitter in value:
            if not re.match(r'^\d{4}-\d{4}-\d{4}-\d{4}$', submitter):
                raise ValueError(f"Submitter must be a valid ORCID, got: {submitter}")
        return value

    @field_validator('overhang_names')
    @classmethod
    def validate_overhang_names(cls, value: Dict[DNASequence, str]) -> Dict[DNASequence, str]:
        for overhang in value.keys():
            if len(overhang) != 4:
                raise ValueError(f"Overhang must be 4 characters long, got: {overhang}")
        return value

    @field_validator('parts')
    @classmethod
    def validate_parts(cls, value: List[Part]) -> List[Part]:
        seen_ids = set()
        for part in value:
            if part.id in seen_ids:
                raise ValueError(f"Duplicate part ID found: {part.id}")
            seen_ids.add(part.id)
        return value

    class Config:
        populate_by_name = True

    def to_edges_graph(self) -> nx.DiGraph:
        graph = nx.DiGraph()
        for part in self.parts:
            graph.add_edge(part.left_overhang, part.right_overhang)
        return graph

    def get_assembly_enzyme(self) -> RestrictionType:
        return parse_restriction_enzymes([self.assembly_enzyme]).format(self.assembly_enzyme)

    def assign_plasmid_to_syntax_part(self, plasmid: Dseqrecord) -> str:
        graph = self.to_edges_graph()
        assembly_enzyme = self.get_assembly_enzyme()
        result = []
        for fragment in plasmid.cut(assembly_enzyme):
            for rc in [True, False]:
                query = fragment.reverse_complement() if rc else fragment
                three_type, three_ovhg = query.seq.three_prime_end()
                five_type, five_ovhg = query.seq.five_prime_end()
                # It must only have 5' overhangs
                if three_type != five_type or five_type != "5'":
                    continue
                # It must not contain the recognition site of the enzyme inside
                # since they are always in the backbone, not the part.
                # We use compsite, because the simple search method requires the
                # cutsite to be there, and not sure how behaviour will be querying
                # the overhangs.
                if assembly_enzyme.compsite.search(str(query.seq)) is not None:
                    continue

                left_node = three_ovhg.upper()
                right_node = reverse_complement(five_ovhg).upper()

                if left_node in graph and right_node in graph and nx.has_path(graph, left_node, right_node):
                    result.append(f"{left_node}-{right_node}")
        return result
