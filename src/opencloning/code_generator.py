from pydna.opencloning_models import CloningStrategy


class CodeGenerator:

    cloning_strategy: CloningStrategy
    import_dict: dict = {
        'pydna.opencloning_models': [
            'CloningStrategy',
        ],
        'pydna.assembly2': [],
        'pydna.dseqrecord': [],
    }

    def __init__(self, cloning_strategy: CloningStrategy):
        self.cloning_strategy = cloning_strategy

    def generate_code(self) -> str:
        """
        Generate python code for the cloning strategy.
        """
        # Find the root sources (no inputs)
        [s for s in self.cloning_strategy.sources if len(s.input) == 0]
